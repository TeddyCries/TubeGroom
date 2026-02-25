import math
import bpy
from mathutils import Vector
from . import utils, geometry

# Local variables for interpolation
tubegroom_data = None
curves_signature = {}

# Face primitive
class Face:
    def __init__(self, subregion_id, center, normal, tangent, bitangent, region_positions):
        self.subregion_id = subregion_id
        self.center = center
        self.normal = normal
        self.region_positions = region_positions
        self.tri_indices = None

# System creation
def create_system(radius=None):
    return {
        'radius': float(radius) if radius is not None else 0.02,
        'curves': [],
        'face_topology': {},
        'region_base': {},
        'curves_map': {}, 
        'global_seed': 42,
    }

# Face building and topology analysis
def build_face(sub_id, positions):
    center = sum(positions, Vector()) / len(positions)
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    n = v1.cross(v2)
    n = n.normalized() if n.length > 1e-9 else Vector((0, 0, 1))
    t = Vector((1, 0, 0))
    b = n.cross(t)
    if b.length < 1e-9:
        b = Vector((0, 1, 0))
    return Face(sub_id, center, n, t, b, positions)
def ensure_root_frame(base, pts, center):
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
def analyze_region_topology(system, region_id):
    region = geometry.TubeGroom.regions.get(region_id)
    if not region:
        return []
    faces = [build_face(sub_id, pos) for sub_id, sub in region.subregions.items() if len(pos := sub.get_positions()) >= 3]
    if not faces:
        return []
    faces.sort(key=lambda x: x.subregion_id)
    root_face = next((f for f in faces if f.subregion_id == 1), faces[0])
    base = system['region_base'].setdefault(region_id, {})
    t, b, n, c = ensure_root_frame(base, root_face.region_positions, root_face.center)
    root_uv = [((p - c).dot(t), (p - c).dot(b)) for p in root_face.region_positions]
    base['root_uv'] = root_uv
    base['n_points'] = len(root_uv)
    tri_idx = utils.tri_face(root_uv)
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

# UV generation and poisson disk sampling
def poisson_disk_hashgrid(poly_uv, tri_indices, radius, seed):
    cs = radius / math.sqrt(2.0)
    minx = maxx = poly_uv[0][0]
    miny = maxy = poly_uv[0][1]
    for p in poly_uv[1:]:
        if p[0] < minx:
            minx = p[0]
        elif p[0] > maxx:
            maxx = p[0]
        if p[1] < miny:
            miny = p[1]
        elif p[1] > maxy:
            maxy = p[1]
    ix0 = int(math.floor(minx / cs))
    ix1 = int(math.ceil(maxx / cs))
    iy0 = int(math.floor(miny / cs))
    iy1 = int(math.ceil(maxy / cs))
    grid = {}
    pts = []
    r_sq = radius * radius
    offsets = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
    for i in range(ix0, ix1):
        for j in range(iy0, iy1):
            a, b = utils.hash2(i, j, seed)
            px = minx + (i - ix0 + a) * cs
            py = miny + (j - iy0 + b) * cs
            if not utils.point_inside_poly((px, py), poly_uv, tri_indices):
                continue
            too_close = False
            for di, dj in offsets:
                if nb := grid.get((i + di, j + dj)):
                    if (px - nb[0])**2 + (py - nb[1])**2 < r_sq:
                        too_close = True
                        break
            if too_close:
                continue
            grid[(i, j)] = (px, py)
            pts.append((px, py))
    return pts
def sample_root_uv(base, radius, seed):
    if not (poly_uv := base.get('root_uv')) or not (tri_indices := base.get('tri_indices')):
        return []
    return poisson_disk_hashgrid(poly_uv, tri_indices, radius, seed)

# Curve generation from region topology
def generate_region_curves(system, region_id):
    region = geometry.TubeGroom.regions.get(region_id)
    if not region:
        return
    base = system['region_base'].setdefault(region_id, {})
    faces = [build_face(sub_id, pos) for sub_id, sub in region.subregions.items() if len(pos := sub.get_positions()) >= 3]
    if not faces:
        return
    faces.sort(key=lambda x: x.subregion_id)
    faces_by_id = {f.subregion_id: f for f in faces}
    sub_ids = sorted(region.subregions.keys())
    root = faces_by_id.get(1, faces[0])
    
    geom_sig = tuple(p.to_tuple(5) for p in root.region_positions)
    if base.get('geom_sig') != geom_sig or base.get('n_points') != len(root.region_positions):
        t, b, n, c = ensure_root_frame(base, root.region_positions, root.center)
        root_uv = [((p - c).dot(t), (p - c).dot(b)) for p in root.region_positions]
        base.update({'root_uv': root_uv, 'n_points': len(root_uv), 'geom_sig': geom_sig})
        tri_idx = utils.tri_face(root_uv)
        if not tri_idx:
            return
        base['tri_indices'] = tri_idx
        base.pop('uv_pts_abs', None)
    else:
        t, b, n, c = base['root_frame']
    
    n_points = base['n_points']
    for f in faces:
        if len(f.region_positions) == n_points:
            f.tri_indices = base['tri_indices']
    system['face_topology'][region_id] = faces
    
    radius = system['radius']
    seed = system['global_seed'] + 1000 * int(region_id)
    
    # Only regenerate UV points if radius changed
    if base.get('last_radius') != radius or 'uv_pts_abs' not in base:
        uv_pts_abs = sample_root_uv(base, radius, seed)
        base.update({'uv_pts_abs': uv_pts_abs, 'last_radius': radius})
    else:
        uv_pts_abs = base['uv_pts_abs']
    
    if not uv_pts_abs:
        return
    
    t_obj = getattr(bpy.context.scene, 'strand_raycast_target', None)
    surface_obj = t_obj if t_obj and t_obj.type == 'MESH' and not utils.tubegroom_object(t_obj) else None
    proj_cache = {}
    
    for point_idx, (u_abs, v_abs) in enumerate(uv_pts_abs):
        curve_points = []
        for sid in sub_ids:
            if not (f := faces_by_id.get(sid)):
                continue
            
            if sid == 1:
                # Root subregion: use original frame + surface projection
                world = c + (t * u_abs) + (b * v_abs)
                if surface_obj:
                    key = (round(u_abs, 5), round(v_abs, 5))
                    if key not in proj_cache:
                        hit, _ = utils.closest_point(surface_obj, world)
                        proj_cache[key] = hit if hit else world
                    world = proj_cache[key]
            else:
                # Non-root subregions: use barycentric interpolation with current positions
                orig_uv = base['root_uv']
                if len(orig_uv) == len(f.region_positions) and len(orig_uv) >= 3:
                    # Find barycentric coordinates of UV point in original triangle space
                    tri_indices = base['tri_indices']
                    world = None
                    
                    # Check each triangle for the point
                    for tri in tri_indices:
                        if len(tri) == 3:
                            i0, i1, i2 = tri
                            if i0 < len(orig_uv) and i1 < len(orig_uv) and i2 < len(orig_uv):
                                p0, p1, p2 = orig_uv[i0], orig_uv[i1], orig_uv[i2]
                                w0, w1, w2, inside = utils.barycentric((u_abs, v_abs), p0, p1, p2)
                                
                                if inside and i0 < len(f.region_positions) and i1 < len(f.region_positions) and i2 < len(f.region_positions):
                                    # Interpolate using current face positions
                                    v0, v1, v2 = f.region_positions[i0], f.region_positions[i1], f.region_positions[i2]
                                    world = v0 * w0 + v1 * w1 + v2 * w2
                                    break
                                
            curve_points.append(world)
        
        if curve_points:
            system['curves'].append({'stream_id': (region_id, point_idx), 'region_id': region_id, 'point_index': point_idx, 'points': curve_points})
def generate_all_curves(system):
    system['curves'] = []
    for rid in sorted(geometry.TubeGroom.regions.keys(), key=lambda r: (int(r) if str(r).isdigit() else float('inf'), str(r))):
        analyze_region_topology(system, rid)
        generate_region_curves(system, rid)
        guide = main_guide_curve(system, rid)
        if guide:
            system['curves'].append(guide)
    system['curves'].sort(key=lambda c: (int(c['region_id']) if str(c['region_id']).isdigit() else float('inf'), str(c['region_id']), c.get('point_index', 0)))
    update_curves_map(system)
    return system['curves']
def generate_region_curves_cached(system, region_id):
    analyze_region_topology(system, region_id)
    system['curves'] = [c for c in system.get('curves', []) if c.get('region_id') != region_id]
    generate_region_curves(system, region_id)
    guide = main_guide_curve(system, region_id)
    if guide:
        system['curves'].append(guide)
    system['curves'].sort(key=lambda c: (int(c.get('region_id')) if str(c.get('region_id')).isdigit() else float('inf'), str(c.get('region_id')), c.get('point_index', 0)))
    update_curves_map(system)

# Curve object building
def update_curves_map(system):
    if not system or 'curves' not in system:
        system['curves_map'] = {}
        return
    from collections import defaultdict
    curves_by_region = defaultdict(list)
    for i, curve in enumerate(system['curves']):
        rid = curve.get('region_id')
        if rid is not None:
            curves_by_region[rid].append(i)
    system['curves_map'] = {rid: {'start_idx': indices[0], 'count': len(indices)} 
                           for rid, indices in curves_by_region.items()}
def main_guide_curve(system, region_id):
    faces = system.get('face_topology', {}).get(region_id, [])
    if not faces:
        return None
    
    # Get surface object for snapping (same as regular curves)
    t_obj = getattr(bpy.context.scene, 'strand_raycast_target', None)
    surface_obj = t_obj if t_obj and t_obj.type == 'MESH' and not utils.tubegroom_object(t_obj) else None
    
    curve_points = []
    for i, face in enumerate(faces):
        barycenter = sum(face.region_positions, Vector()) / len(face.region_positions)
        # Only snap the first point to surface (like regular curves do with root subregion)
        if i == 0 and surface_obj:
            hit, _ = utils.closest_point(surface_obj, barycenter)
            if hit:
                barycenter = hit
        curve_points.append(barycenter)
    max_idx = max((c.get('point_index', -1) for c in system.get('curves', []) if c.get('region_id') == region_id), default=-1)
    return {'stream_id': (region_id, 'guide'), 'region_id': region_id, 'point_index': max_idx + 1, 'points': curve_points, 'is_guide': True}
def build_curves_object(base_name, system):
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
    global curves_signature
    sig_map = curves_signature
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
    update_curves_map(system)
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
    curves_signature[name_obj] = curr_sig
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
def generate_interpolation():
    global tubegroom_data
    if not geometry.TubeGroom.regions:
        system = create_system()
        tubegroom_data = system
        return system
    
    current_radius = float(getattr(bpy.context.scene, 'strand_poisson_radius', 0.02))
    
    if isinstance(tubegroom_data, dict) and 'curves' in tubegroom_data:
        tubegroom_data['radius'] = current_radius
        old_radius = tubegroom_data.get('last_global_radius')
        if old_radius is None or abs(old_radius - current_radius) > 1e-6:
            tubegroom_data['last_global_radius'] = current_radius
            for base in tubegroom_data.get('region_base', {}).values():
                base.pop('uv_pts_abs', None)
                base.pop('last_radius', None)
        generate_all_curves(tubegroom_data)
        return tubegroom_data
    
    # Create new system
    system = create_system(radius=current_radius)
    system['last_global_radius'] = current_radius
    generate_all_curves(system)
    tubegroom_data = system
    return system

# Interpolation management functions
def update_interpolation(context, region_id=None, update_topology=False):
    if not getattr(context.scene, 'strand_curves_enabled', False):
        return
    
    global tubegroom_data
    base_obj = geometry.get_tg_obj(context, allow_create=False)
    if not base_obj:
        return
    
    base_name = utils.get_base_name(base_obj)
    if not base_name:
        return
    
    if not geometry.TubeGroom.regions:
        tubegroom_data = None
        system = create_system()
        build_curves_object(base_name, system)
        return
    
    if not getattr(context.scene, 'strand_interpolation_enabled', False):
        return
    
    system = tubegroom_data or generate_interpolation()
    if not isinstance(system, dict):
        system = generate_interpolation()
    
    current_radius = float(getattr(context.scene, 'strand_poisson_radius', 0.02))
    if system.get('radius') != current_radius:
        system['radius'] = current_radius
        for base in system.get('region_base', {}).values():
            base.pop('uv_pts_abs', None)
            base.pop('last_radius', None)
    
    if region_id is None or update_topology:
        existing_region_ids = set(geometry.TubeGroom.regions.keys())
        for rid in list(system.get('face_topology', {}).keys()):
            if rid not in existing_region_ids:
                system['face_topology'].pop(rid, None)
        for rid in list(system.get('region_base', {}).keys()):
            if rid not in existing_region_ids:
                system['region_base'].pop(rid, None)
        generate_all_curves(system)
    else:
        generate_region_curves_cached(system, region_id)
    
    tubegroom_data = system
    curves_obj = build_curves_object(base_name, system)
    if curves_obj:
        live_enabled = getattr(context.scene, 'strand_interpolation_enabled', False)
        curves_obj.show_in_front = live_enabled
def rebuild_regions(obj):
    if not obj or obj.type != 'MESH':
        return False
    mesh = obj.data
    if not all(attr in mesh.attributes for attr in ['region_id', 'subregion_id']):
        return False
    geometry.TubeGroom.regions.clear()
    geometry.TubeGroom.next_region_id = 1
    from collections import defaultdict
    verts_by_subregion = defaultdict(list)
    rid_attr = mesh.attributes['region_id'].data
    sid_attr = mesh.attributes['subregion_id'].data
    for i in range(len(mesh.vertices)):
        rid = rid_attr[i].value
        sid = sid_attr[i].value
        if rid > 0 and sid > 0:
            verts_by_subregion[(rid, sid)].append(i)
    sorted_points_by_subregion = {}
    num_verts = len(mesh.vertices)
    adj_list = [[] for _ in range(num_verts)]
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
            matrix = obj.matrix_world
            sorted_points_by_subregion[(rid, sid)] = [
                matrix @ mesh.vertices[v_idx].co
                for v_idx in ordered_indices
            ]
    temp_regions = {}
    for (rid, sid), points in sorted_points_by_subregion.items():
        if rid not in temp_regions:
            temp_regions[rid] = {}
        temp_regions[rid][sid] = points
    max_rid = 0
    for rid, region_data in temp_regions.items():
        if rid > max_rid:
            max_rid = rid
        region = geometry.Region(rid, rid)
        max_sid = 0
        for sid, points in sorted(region_data.items()):
            if sid > max_sid:
                max_sid = sid
            subregion = geometry.Subregion(rid, sid)
            for pos in points:
                subregion.add_point(pos)
            region.subregions[sid] = subregion
        region.next_subregion_id = max_sid + 1
        geometry.TubeGroom.regions[rid] = region
    geometry.TubeGroom.next_region_id = max_rid + 1
    # Color index managed internally
    return len(geometry.TubeGroom.regions) > 0
