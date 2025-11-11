import math
import bpy
from mathutils import Vector
from . import core as shared_data, utils
from .core import Face

EPS = 1e-9

# Face primitive
def tubegroom_default_radius():
    val = getattr(shared_data, 'poisson_radius', None)
    if val is not None:
        return float(val)
    return float(getattr(bpy.context.scene, 'strand_poisson_radius', 0.02))
def tubegroom_system_create(radius=None):
    return {
        'radius': float(radius) if radius is not None else tubegroom_default_radius(),
        'curves': [],
        'face_topology': {},
        'region_base': {},
        'curves_map': {}, # region_id -> {'start_idx': int, 'count': int}
        'global_seed': 42,
    }
def tubegroom_build_face(sub_id, positions):
    center = sum(positions, Vector()) / len(positions)
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    n = v1.cross(v2)
    n = n.normalized() if n.length > EPS else Vector((0, 0, 1))
    t = Vector((1, 0, 0))
    b = n.cross(t)
    if b.length < EPS:
        b = Vector((0, 1, 0))
    return Face(sub_id, center, n, t, b, positions)
def tubegroom_ensure_root_frame(base, pts, center):
    if 'root_frame' in base:
        return base['root_frame']
    acc = Vector()
    m = len(pts)
    for i in range(m):
        acc += (pts[i] - pts[i - 1]).cross(pts[(i + 1) % m] - pts[i])
    n = acc.normalized() if acc.length > 0 else Vector((0, 0, 1))
    up = Vector((0, 0, 1)) if abs(n.z) < 0.99 else Vector((0, 1, 0))
    t = (up.cross(n)).normalized() if up.cross(n).length else Vector((1, 0, 0))
    b = (n.cross(t)).normalized()
    base['root_frame'] = (t, b, n, center.copy())
    return base['root_frame']
def tubegroom_analyze_region_topology(system, region_id):
    region = shared_data.regions.get(region_id)
    if not region:
        return []
    faces = [tubegroom_build_face(sub_id, pos) for sub_id, sub in region.subregions.items() if len(pos := sub.get_positions()) >= 3]
    if not faces:
        return []
    faces.sort(key=lambda x: x.subregion_id)
    root_face = next((f for f in faces if f.subregion_id == 1), faces[0])
    base = system['region_base'].setdefault(region_id, {})
    t, b, n, c = tubegroom_ensure_root_frame(base, root_face.region_positions, root_face.center)
    root_uv = [((p - c).dot(t), (p - c).dot(b)) for p in root_face.region_positions]
    base['root_uv'] = root_uv
    base['n_points'] = len(root_uv)
    tri_idx = utils.triangulate_face(root_uv)
    if not tri_idx:
        return []
    base['tri_indices'] = tri_idx
    n_points = base['n_points']
    for f in faces:
        if len(f.region_positions) == n_points:
            f.reference_vertices2d = root_uv
            f.tri_indices = tri_idx
    system['face_topology'][region_id] = faces
    return faces
def tubegroom_poisson_disk_hashgrid(poly_uv, tri_indices, radius, seed):
    cs = radius / math.sqrt(2.0)
    xs = [p[0] for p in poly_uv]
    ys = [p[1] for p in poly_uv]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    ix0 = int(math.floor(minx / cs))
    ix1 = int(math.ceil(maxx / cs))
    iy0 = int(math.floor(miny / cs))
    iy1 = int(math.ceil(maxy / cs))
    grid = {}
    pts = []
    r_sq = radius * radius
    for i in range(ix0, ix1):
        for j in range(iy0, iy1):
            a, b = utils.hash2(i, j, seed)
            px = minx + (i - ix0 + a) * cs
            py = miny + (j - iy0 + b) * cs
            if not utils.point_in_triangulated_poly((px, py), poly_uv, tri_indices):
                continue
            if any((px - nb[0])**2 + (py - nb[1])**2 < r_sq for di in (-1, 0, 1) for dj in (-1, 0, 1) if (nb := grid.get((i + di, j + dj)))):
                continue
            grid[(i, j)] = (px, py)
            pts.append((px, py))
    return pts
def tubegroom_uv_sample_root(base, radius, seed):
    if not (poly_uv := base.get('root_uv')) or not (tri_indices := base.get('tri_indices')):
        return []
    return tubegroom_poisson_disk_hashgrid(poly_uv, tri_indices, radius, seed)
def tubegroom_generate_curves_for_region(system, region_id):
    region = shared_data.regions.get(region_id)
    if not region:
        return
    base = system['region_base'].setdefault(region_id, {})
    faces = [tubegroom_build_face(sub_id, pos) for sub_id, sub in region.subregions.items() if len(pos := sub.get_positions()) >= 3]
    if not faces:
        return
    faces.sort(key=lambda x: x.subregion_id)
    faces_by_id = {f.subregion_id: f for f in faces}
    sub_ids = sorted(region.subregions.keys())
    root = faces_by_id.get(1, faces[0])
    t, b, n, c = tubegroom_ensure_root_frame(base, root.region_positions, root.center)
    geom_sig = tuple(p.to_tuple(5) for p in root.region_positions)
    if base.get('geom_sig') != geom_sig or base.get('n_points') != len(root.region_positions):
        root_uv = [((p - c).dot(t), (p - c).dot(b)) for p in root.region_positions]
        base.update({'root_uv': root_uv, 'n_points': len(root_uv), 'geom_sig': geom_sig})
        tri_idx = utils.triangulate_face(root_uv)
        if not tri_idx:
            return
        base['tri_indices'] = tri_idx
        base.pop('uv_pts_abs', None)
    n_points = base['n_points']
    for f in faces:
        if len(f.region_positions) == n_points:
            f.reference_vertices2d = base['root_uv']
            f.tri_indices = base['tri_indices']
    system['face_topology'][region_id] = faces
    radius = system['radius']
    seed = system['global_seed'] + 1000 * int(region_id)
    uv_sig = (len(base['root_uv']), tuple((round(u, 6), round(v, 6)) for (u, v) in base['root_uv']))
    if base.get('uv_sig') == uv_sig and base.get('last_radius') == radius and 'uv_pts_abs' in base:
        uv_pts_abs = base['uv_pts_abs']
    else:
        uv_pts_abs = tubegroom_uv_sample_root(base, radius, seed)
        base.update({'uv_pts_abs': uv_pts_abs, 'uv_sig': uv_sig, 'last_radius': radius})
    if not uv_pts_abs:
        return
    t_obj = getattr(bpy.context.scene, 'strand_raycast_target', None)
    surface_obj = t_obj if t_obj and t_obj.type == 'MESH' and not utils.is_tubegroom_object(t_obj) else None
    proj_cache = {}
    for point_idx, (u_abs, v_abs) in enumerate(uv_pts_abs):
        curve_points = []
        for sid in sub_ids:
            if not (f := faces_by_id.get(sid)):
                continue
            world = f.uv_to_world_abs(u_abs, v_abs)
            if sid == 1 and surface_obj:
                key = (round(u_abs, 5), round(v_abs, 5))
                if key not in proj_cache:
                    hit, _ = utils._closest_point_on_surface_world(surface_obj, world)
                    proj_cache[key] = hit if hit else world
                world = proj_cache[key]
            curve_points.append(world)
        if curve_points:
            system['curves'].append({'stream_id': (region_id, point_idx), 'region_id': region_id, 'point_index': point_idx, 'points': curve_points})
def tubegroom_generate_all_curves(system):
    system['curves'] = []
    for rid in sorted(shared_data.regions.keys(), key=lambda r: (int(r) if str(r).isdigit() else float('inf'), str(r))):
        tubegroom_analyze_region_topology(system, rid)
        tubegroom_generate_curves_for_region(system, rid)
        guide = main_guide_curve(system, rid)
        if guide:
            system['curves'].append(guide)
    system['curves'].sort(key=lambda c: (int(c['region_id']) if str(c['region_id']).isdigit() else float('inf'), str(c['region_id']), c.get('point_index', 0)))
    _update_curves_map(system)
    return system['curves']
def tubegroom_generate_curves_for_region_cached(system, region_id):
    tubegroom_analyze_region_topology(system, region_id)
    system['curves'] = [c for c in system.get('curves', []) if c.get('region_id') != region_id]
    tubegroom_generate_curves_for_region(system, region_id)
    guide = main_guide_curve(system, region_id)
    if guide:
        system['curves'].append(guide)
    system['curves'].sort(key=lambda c: (int(c.get('region_id')) if str(c.get('region_id')).isdigit() else float('inf'), str(c.get('region_id')), c.get('point_index', 0)))
    _update_curves_map(system)
def _update_curves_map(system):
    if not system or 'curves' not in system:
        system['curves_map'] = {}
        return
    curves_by_region = {}
    for i, curve in enumerate(system['curves']):
        rid = curve.get('region_id')
        if rid is not None:
            curves_by_region.setdefault(rid, []).append(i)
    system['curves_map'] = {rid: {'start_idx': min(indices), 'count': len(indices)} for rid, indices in curves_by_region.items() if indices}
def main_guide_curve(system, region_id):
    faces = system.get('face_topology', {}).get(region_id, [])
    if not faces:
        return None
    curve_points = []
    for face in faces:
        barycenter = sum(face.region_positions, Vector()) / len(face.region_positions)
        curve_points.append(barycenter)
    max_idx = max((c.get('point_index', -1) for c in system.get('curves', []) if c.get('region_id') == region_id), default=-1)
    return {'stream_id': (region_id, 'guide'), 'region_id': region_id, 'point_index': max_idx + 1, 'points': curve_points, 'is_guide': True}
def build_curves_object_from_system(base_name, system):
    name_obj = f"CRV_{base_name}"
    streams = []
    curves = system.get('curves', []) if isinstance(system, dict) else []
    if curves:
        streams = [c for c in curves if c.get('points') and len(c['points']) >= 2]
        def _k(c):
            rid = c.get('region_id')
            rk = (0, int(rid)) if str(rid).isdigit() else (1, str(rid))
            return (rk, int(c.get('point_index', 0)))
        streams.sort(key=_k)
    obj = bpy.data.objects.get(name_obj)
    old_data = obj.data if obj else None
    new_sizes = [len(c['points']) for c in streams]
    curr_sig = tuple(sorted([(s.get('region_id'), len(s['points'])) for s in streams]))
    sig_map = getattr(shared_data, 'curves_signature', {}) or {}
    prev_sig = sig_map.get(name_obj)
    topology_changed = (prev_sig != curr_sig)
    if not topology_changed and old_data and hasattr(old_data, 'curves') and len(old_data.curves) == len(new_sizes):
        old_sizes = [len(c.points) for c in old_data.curves]
        if old_sizes == new_sizes:
            coords = []
            for s in streams:
                for p in s['points']:
                    coords.extend([p.x, p.y, p.z])
            if 'position' in old_data.attributes:
                old_data.attributes["position"].data.foreach_set("vector", coords)
            if hasattr(old_data, 'update_tag'):
                old_data.update_tag()
            return obj
    _update_curves_map(system)
    curves_db = bpy.data.hair_curves
    new_data = curves_db.new(name=name_obj + "_data_tmp")
    if streams:
        new_data.add_curves(new_sizes)
        coords = []
        for s in streams:
            for p in s['points']:
                coords.extend([p.x, p.y, p.z])
        new_data.attributes["position"].data.foreach_set("vector", coords)
        if 'guide_curve_index' not in new_data.attributes:
            new_data.attributes.new('guide_curve_index', 'INT', 'CURVE')
        guide_indices = []
        for i, s in enumerate(streams):
            rid = s.get('region_id')
            guide_idx = next((j for j, c in enumerate(streams) if c.get('region_id') == rid and c.get('is_guide')), i)
            guide_indices.append(guide_idx)
        new_data.attributes['guide_curve_index'].data.foreach_set('value', guide_indices)
    if old_data:
        for mat in old_data.materials:
            if mat:
                new_data.materials.append(mat)
    if not hasattr(shared_data, 'curves_signature'):
        shared_data.curves_signature = {}
    shared_data.curves_signature[name_obj] = curr_sig
    if obj is None:
        obj = bpy.data.objects.new(name_obj, new_data)
        if hasattr(bpy.context.scene, 'collection'):
            bpy.context.scene.collection.objects.link(obj)
        elif hasattr(bpy.context, 'collection'):
            bpy.context.collection.objects.link(obj)
    else:
        obj.data = new_data
    if old_data and old_data.users == 0 and hasattr(curves_db, 'remove'):
        curves_db.remove(old_data)
    if hasattr(new_data, 'update_tag'):
        new_data.update_tag()
    if hasattr(obj, 'hide_set'):
        obj.hide_set(False)
    if hasattr(obj, 'hide_viewport'):
        obj.hide_viewport = False
    if hasattr(obj, 'hide_render'):
        obj.hide_render = False
        
    return obj
def generate_tubegroom_interpolation():
    if not shared_data.regions:
        system = tubegroom_system_create()
        shared_data.tubegroom_data = system
        return system
    r = getattr(shared_data, 'poisson_radius', None)
    if r is None and hasattr(bpy.context.scene, 'strand_poisson_radius'):
        r = float(getattr(bpy.context.scene, 'strand_poisson_radius'))
    if r is None:
        r = 0.02
    sys_prev = getattr(shared_data, 'tubegroom_data', None)
    if isinstance(sys_prev, dict) and 'curves' in sys_prev:
        sys_prev['radius'] = float(r)
        tubegroom_generate_all_curves(sys_prev)
        return sys_prev
    system = tubegroom_system_create(radius=float(r))
    tubegroom_generate_all_curves(system)
    shared_data.tubegroom_data = system
    return system

# Interpolation management functions
def update_tubegroom_interpolation(context, region_id=None, update_topology=False):
    if not getattr(context.scene, 'tubegroom_curves_enabled', False):
        return
    if not getattr(context.scene, 'strand_interpolation_enabled', False):
        return
    if not shared_data.regions:
        shared_data.tubegroom_data = None
        from . import geometry
        base_obj = geometry.get_or_create_tubegroom_object(context, allow_create=False)
        if base_obj:
            base_name = utils.get_base_name_from_obj(base_obj)
            if base_name:
                system = tubegroom_system_create()
                build_curves_object_from_system(base_name, system)
        if hasattr(shared_data, 'curves_signature'):
            shared_data.curves_signature.pop(f"CRV_{utils.get_base_name_from_obj(base_obj)}" if base_obj else None, None)
        return
    system = shared_data.tubegroom_data or generate_tubegroom_interpolation()
    if not isinstance(system, dict):
        system = generate_tubegroom_interpolation()
    
    if region_id is None or update_topology:
        existing_region_ids = set(shared_data.regions.keys())
        for rid in list(system.get('face_topology', {}).keys()):
            if rid not in existing_region_ids:
                system['face_topology'].pop(rid, None)
        for rid in list(system.get('region_base', {}).keys()):
            if rid not in existing_region_ids:
                system['region_base'].pop(rid, None)
        tubegroom_generate_all_curves(system)
    else:
        tubegroom_generate_curves_for_region_cached(system, region_id)
    shared_data.tubegroom_data = system
    from . import geometry
    
    base_obj = geometry.get_or_create_tubegroom_object(context, allow_create=False)
    if not base_obj:
        return
    
    base_name = utils.get_base_name_from_obj(base_obj)
    if not base_name:
        return
    
    curves_obj = build_curves_object_from_system(base_name, system)
    if curves_obj:
        curves_obj.show_in_front = True
def rebuild_regions_from_merged_mesh(obj):
    if not obj or obj.type != 'MESH':
        return False
    scenes = getattr(obj, 'users_scene', None)
    scene = scenes[0] if scenes else None
    shared_data.set_history_owner(obj, scene)
    mesh = obj.data
    if not all(attr in mesh.attributes for attr in ['region_id', 'subregion_id']):
        return False
    shared_data.reset_all_data(clear_history=False)
    # Step 1: Group vertex indices by (rid, sid)
    verts_by_subregion = {}
    rid_attr = mesh.attributes['region_id'].data
    sid_attr = mesh.attributes['subregion_id'].data
    for i, v in enumerate(mesh.vertices):
        rid = rid_attr[i].value
        sid = sid_attr[i].value
        if rid <= 0 or sid <= 0:
            continue
        key = (rid, sid)
        if key not in verts_by_subregion:
            verts_by_subregion[key] = []
        verts_by_subregion[key].append(i)
    # Step 2: For each subregion, sort vertices by walking the edges
    sorted_points_by_subregion = {}
    # Create an adjacency list for the entire mesh once
    adj_list = [[] for _ in range(len(mesh.vertices))]
    for edge in mesh.edges:
        v1, v2 = edge.vertices
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    for (rid, sid), vert_indices in verts_by_subregion.items():
        if len(vert_indices) < 3:
            continue
        vert_set = set(vert_indices)
        # Find a starting vertex (one with exactly two neighbors within the subregion loop)
        start_node = -1
        for v_idx in vert_indices:
            # Count how many neighbors of this vertex are also in the same subregion
            subregion_neighbors = [n for n in adj_list[v_idx] if n in vert_set]
            if len(subregion_neighbors) == 2:
                start_node = v_idx
                break
        if start_node == -1:
            continue
        # Walk the loop using the adjacency list
        ordered_indices = []
        curr, prev = start_node, -1
        while len(ordered_indices) < len(vert_indices):
            ordered_indices.append(curr)
            neighbors = [n for n in adj_list[curr] if n in vert_set]
            next_node = next((n for n in neighbors if n != prev), None)
            if next_node is None:
                break
            prev, curr = curr, next_node
        if len(ordered_indices) == len(vert_indices):
            sorted_points_by_subregion[(rid, sid)] = [
                obj.matrix_world @ mesh.vertices[v_idx].co
                for v_idx in ordered_indices
            ]
    # Step 3: Rebuild the core.Region and core.Subregion objects
    temp_regions = {}
    for (rid, sid), points in sorted_points_by_subregion.items():
        if rid not in temp_regions:
            temp_regions[rid] = {}
        temp_regions[rid][sid] = points
    max_rid = 0
    for rid, region_data in temp_regions.items():
        if rid > max_rid:
            max_rid = rid
        region = shared_data.Region(rid, rid - 1)
        max_sid = 0
        for sid, points in sorted(region_data.items()):
            if sid > max_sid:
                max_sid = sid
            subregion = shared_data.Subregion(rid, sid)
            for pos in points:
                subregion.add_point(pos)
            region.subregions[sid] = subregion
        region.next_subregion_id = max_sid + 1
        shared_data.regions[rid] = region
    shared_data.next_region_id = max_rid + 1
    shared_data.current_region_color_index = max_rid
    return len(shared_data.regions) > 0
