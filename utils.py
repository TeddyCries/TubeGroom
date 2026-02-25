import bpy
import colorsys
from mathutils import Vector
from bpy_extras import view3d_utils
from collections import deque
from . import geometry

# Color utilities
def color_id(idx):
    """Generates a visually distinct color based on an index using the golden ratio."""
    golden_ratio_conjugate = 0.61803398875
    h = (idx * golden_ratio_conjugate) % 1.0
    r, g, b = colorsys.hsv_to_rgb(h, 1.0, 1.0)
    return (r, g, b, 0.5)
def regions_color(positions, colors):
    """Group positions by their corresponding RGBA color tuples."""
    if not positions or not colors:
        return {}
    groups = {}
    for pos, col in zip(positions, colors):
        key = tuple(col) if isinstance(col, (list, tuple)) else col
        groups.setdefault(key, []).append(pos)
    return groups

# Object naming utilities
def tubegroom_object(obj):
    """Checks if an object is a TubeGroom object by its name prefix."""
    return obj and obj.name.startswith("GEO_TubeGroom_")
def get_base_name(obj):
    """Extracts the base name (e.g., "TubeGroom_Cube") from a GEO object."""
    if not tubegroom_object(obj):
        return None
    return obj.name[14:]

# Math and geometry utilities
def surface_normal(obj, point):
    """Calculates the world-space surface normal of an object at a given world-space point."""
    _, normal = closest_point(obj, point)
    return normal if normal is not None else Vector((0, 0, 1))
def project_2d(point_3d, context):
    """Projects a 3D world-space point to 2D screen coordinates."""
    return view3d_utils.location_3d_to_region_2d(context.region, context.region_data, point_3d)
def dist2_seg3(p, a, b):
    """Squared distance from point p to 3D segment a-b."""
    ab = b - a
    l2 = ab.length_squared
    if l2 == 0.0:
        return (p - a).length_squared
    t = max(0, min(1, (p - a).dot(ab) / l2))
    return (p - (a + t * ab)).length_squared
def dist2_seg2(point, line_start, line_end):
    """Calculates the squared distance from a point to a 2D line segment."""
    px, py = point
    x1, y1 = line_start
    x2, y2 = line_end

    c, d = x2 - x1, y2 - y1
    len_sq = c * c + d * d
    if len_sq == 0:
        return (px - x1) ** 2 + (py - y1) ** 2

    t = max(0, min(1, ((px - x1) * c + (py - y1) * d) / len_sq))
    xx, yy = x1 + t * c, y1 + t * d
    return (px - xx) ** 2 + (py - yy) ** 2
def barycentric(p, a, b, c):
    """Calculates barycentric coordinates of a point p inside a triangle (a, b, c)."""
    px, py = p
    ax, ay = a
    bx, by = b
    cx, cy = c

    v0x, v0y = (cx - ax), (cy - ay)
    v1x, v1y = (bx - ax), (by - ay)
    v2x, v2y = (px - ax), (py - ay)

    dot00 = v0x*v0x + v0y*v0y
    dot01 = v0x*v1x + v0y*v1y
    dot02 = v0x*v2x + v0y*v2y
    dot11 = v1x*v1x + v1y*v1y
    dot12 = v1x*v2x + v1y*v2y

    denom = dot00 * dot11 - dot01 * dot01

    if abs(denom) < 1e-20:
        return (1/3, 1/3, 1/3, False)
    inv = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    w = 1.0 - u - v
    ok = (u >= -1e-6) and (v >= -1e-6) and (w >= -1e-6)
    return (w, v, u, ok)
def tri_face(points, tol=1e-9):
    """Triangulate a polygon loop using ear clipping."""
    n = len(points)
    if n < 3:
        return []
    uv = [(float(p[0]), float(p[1])) for p in points]
    idx = list(range(n))
    def area2(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    a = 0.0
    for i in range(n):
        x1, y1 = uv[i]
        x2, y2 = uv[(i+1)%n]
        a += x1*y2 - x2*y1
    if a < 0:
        uv.reverse()
        idx.reverse()
    def point_in_tri(p, a, b, c):
        s1 = area2(a,b,p)
        s2 = area2(b,c,p)
        s3 = area2(c,a,p)

        has_neg = (s1 < -tol) or (s2 < -tol) or (s3 < -tol)
        has_pos = (s1 > tol) or (s2 > tol) or (s3 > tol)
        return not (has_neg and has_pos)
    res = []
    v = list(range(len(uv)))
    guard = 0
    while len(v) > 3 and guard < len(points) * len(points):
        guard += 1
        m = len(v)
        for k in range(m):
            i0, i1, i2 = v[(k-1)%m], v[k], v[(k+1)%m]
            a, b, c = uv[i0], uv[i1], uv[i2]
            if area2(a, b, c) <= tol:
                continue
            ok = True
            for j in v:
                if j in (i0, i1, i2):
                    continue
                if point_in_tri(uv[j], a, b, c):
                    ok = False
                    break
            if not ok:
                continue
            # Clip ear
            res.append([idx[i0], idx[i1], idx[i2]])
            del v[k]
            break
    if len(v) == 3:
        res.append([idx[v[0]], idx[v[1]], idx[v[2]]])
    return res
def tip_normal(positions, surface_obj):
    """Return the average normal of a face loop, oriented away from the surface."""
    if len(positions) < 3:
        return Vector((0, 0, 1))

    normal = (positions[1] - positions[0]).cross(positions[2] - positions[0])
    if normal.length == 0:
        return Vector((0, 0, 1))
    
    normal.normalize()
    if surface_obj and surface_obj.type == 'MESH':
        surf_n = surface_normal(surface_obj, sum(positions, Vector()) / len(positions))
        if surf_n and surf_n.length > 0 and normal.dot(surf_n) < 0:
            normal.negate()
    return normal

# Mesh and topology utilities
def opposite_edge(face, edge):
    """In a quad face, return the edge opposite to the given edge."""
    for ed in face.edges:
        if ed is edge:
            continue

        va = ed.verts[0]
        vb = ed.verts[1]

        if (va not in edge.verts) and (vb not in edge.verts):
            return ed
    return None
def find_islands(mesh_obj):
    """Return list of sets of vertex indices representing connected components."""
    mesh = mesh_obj.data
    n = len(mesh.vertices)
    if n == 0:
        return []
    adj = [[] for _ in range(n)]
    for e in mesh.edges:
        i, j = e.vertices[0], e.vertices[1]
        adj[i].append(j)
        adj[j].append(i)
    visited = [False] * n
    islands = []
    for i in range(n):
        if visited[i]:
            continue
        # Breadth-first search
        comp = set()
        queue = [i]
        visited[i] = True
        while queue:
            u = queue.pop(0)
            comp.add(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    queue.append(v)
        islands.append(comp)
    return islands
def boundary_loops(bm, allowed_idx=None):
    """Finds all distinct boundary loops in a BMesh, returning sets of vertex indices."""
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()

    allowed = None if allowed_idx is None else set(allowed_idx)
    boundary_edges = [e for e in bm.edges if e.is_boundary and (allowed is None or (e.verts[0].index in allowed and e.verts[1].index in allowed))]
    if not boundary_edges:
        return []
    adj = {}
    for e in boundary_edges:
        i0, i1 = e.verts[0].index, e.verts[1].index
        adj.setdefault(i0, set()).add(i1)
        adj.setdefault(i1, set()).add(i0)
    unvisited = set(adj.keys())
    loops = []
    while unvisited:
        start = unvisited.pop()
        stack = [start]
        comp = {start}
        while stack:
            i = stack.pop()
            for j in adj.get(i, ()): 
                if j in unvisited:
                    unvisited.remove(j)
                    comp.add(j)
                    stack.append(j)
        loops.append(comp)
    return loops
def closest_point(obj, world_point):
    """Return (hit_world, world_normal)"""
    if not obj or obj.type != 'MESH':
        return None, None
    depsgraph = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(depsgraph)
    mat = obj_eval.matrix_world
    inv = mat.inverted()
    ok, loc, normal, *_ = obj_eval.closest_point_on_mesh(inv @ world_point)
    if not ok:
        return None, None
    hit_world = mat @ loc
    world_normal = (mat.to_3x3() @ normal).normalized()
    return hit_world, world_normal
def subdiv_edges(positions, subdivisions=2):
    """Interpolates points along the edges of a polygon loop."""
    if len(positions) < 2 or subdivisions <= 1:
        return positions
    subdivided = []
    for i in range(len(positions)):
        current = positions[i]
        next_pos = positions[(i + 1) % len(positions)]
        subdivided.append(current)
        for j in range(1, subdivisions):
            t = j / subdivisions
            subdivided.append(current.lerp(next_pos, t))
    return subdivided
def snap_points(obj, positions, offset=0.001, subdivisions=2):
    """Projects a loop of points onto a surface object with an offset."""
    if not obj or len(positions) < 3:
        return positions
    subdivided_positions = subdiv_edges(positions, subdivisions)
    shrink_points = []
    for pos in subdivided_positions:
        if not obj.data or not obj.data.vertices:
            shrink_points.append(pos)
            continue
        hit_world, world_normal = closest_point(obj, pos)
        if hit_world is not None and world_normal is not None:
            shrink_points.append(hit_world + (world_normal * float(offset)))
        else:
            shrink_points.append(pos)
    return shrink_points
def order_ring_from_indices(bm, idx_set):
    """Orders a set of vertex indices into a single loop by walking edges within the set."""
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    idx_set = set(idx_set)
    if not idx_set:
        return []
    adj = {i: [] for i in idx_set}
    for e in bm.edges:
        a, b = e.verts[0].index, e.verts[1].index
        if a in idx_set and b in idx_set:
            adj[a].append(b)
            adj[b].append(a)
    start = next((i for i in idx_set if len(adj[i]) == 2), next(iter(idx_set), None))
    if start is None:
        return []
    ordered = [start]
    prev, curr = None, start
    for _ in range(len(idx_set)):
        neighbors = adj.get(curr, [])
        nxt = next((nb for nb in neighbors if nb != prev), None)
        if nxt is None or (nxt == ordered[0] and len(set(ordered)) >= len(idx_set)):
            break
        ordered.append(nxt)
        prev, curr = curr, nxt
    return ordered if len(set(ordered)) >= len(idx_set) else []

# Sampling and polygon utilities
def hash2(i, j, seed):
    """Simple 2D integer hash function returning two pseudo-random floats [0, 1]."""
    x = (i * 73856093) ^ (j * 19349663) ^ (seed * 83492791)
    x ^= (x >> 13) & 0xFFFFFFFF
    x = (x * 1274126177) & 0xFFFFFFFF
    x ^= (x >> 16)
    a = (x & 0xFFFF) / 65535.0
    b = ((x >> 16) & 0xFFFF) / 65535.0
    return a, b
def point_inside_poly(p, poly_uv, tri_indices):
    """Checks if a 2D point is inside a pre-triangulated polygon."""
    for (i, j, k) in tri_indices:

        a = poly_uv[i]
        b = poly_uv[j]
        c = poly_uv[k]

        _, _, _, inside = barycentric(p, a, b, c)
        if inside:
            return True
    return False

# Bend interpolation utilities
def get_bend_and_segs(context):
    """Return bend factor and view segments, clamped to valid ranges."""
    scene = context.scene
    bend = max(0.0, min(1.0, getattr(scene, 'strand_bend_factor', 0.0)))
    segs = max(0, getattr(scene, 'strand_view_segments', 0))
    return bend, segs
def interp_rings(real_rings, bend, segs):
    """Create intermediate rings between real rings using Hermite interpolation."""
    if not real_rings or segs <= 0:
        return real_rings
    
    last_index = len(real_rings) - 1
    if last_index == 0:
        return real_rings
    
    inv_bend = 1.0 - bend
    rings = [real_rings[0]]
    segs_plus_1 = segs + 1
    inv_segs = 1.0 / segs_plus_1
    
    for i in range(last_index):
        ring_a_pts, meta_a = real_rings[i]
        ring_b_pts, _ = real_rings[i + 1]
        n = len(ring_a_pts)
        
        if n < 3 or n != len(ring_b_pts):
            continue
        
        prev_ring_pts = real_rings[i - 1][0] if i > 0 and len(real_rings[i - 1][0]) == n else ring_a_pts
        next_ring_pts = real_rings[i + 2][0] if i + 2 <= last_index and len(real_rings[i + 2][0]) == n else ring_b_pts
        
        tangents_a = [(pb - pm) * 0.5 for pb, pm in zip(ring_b_pts, prev_ring_pts)]
        tangents_b = [(pn - pa) * 0.5 for pa, pn in zip(ring_a_pts, next_ring_pts)]
        
        for k in range(1, segs_plus_1):
            t = k * inv_segs
            t2 = t * t
            t3 = t2 * t
            h00 = 2.0 * t3 - 3.0 * t2 + 1.0
            h10 = t3 - 2.0 * t2 + t
            h01 = -2.0 * t3 + 3.0 * t2
            h11 = t3 - t2
            
            h00_bend = h00 * bend
            h10_bend = h10 * bend
            h01_bend = h01 * bend
            h11_bend = h11 * bend
            t_inv = t * inv_bend
            omt_inv = (1.0 - t) * inv_bend
            
            loop = [pa * (h00_bend + omt_inv) + ta * h10_bend + pb * (h01_bend + t_inv) + tb * h11_bend
                    for pa, pb, ta, tb in zip(ring_a_pts, ring_b_pts, tangents_a, tangents_b)]
            
            rings.append((loop, (meta_a[0], 0)))
        
        if i < last_index - 1:
            rings.append(real_rings[i + 1])
    
    rings.append(real_rings[-1])
    return rings

# Mouse interaction utilities
def get_surface_obj(context):
    """Return the surface object selected in the UI picker, if valid."""
    target = getattr(context.scene, 'strand_raycast_target', None)
    return target if target and target.type == 'MESH' and not tubegroom_object(target) else None
def ray_surface(context, mouse_2d, obj, avoid_tg=True):
    """Perform raycast from mouse to surface object."""
    if not obj:
        return None
    if avoid_tg:
        # Check if ray hits TubeGroom object first
        tg_hit, _ = ray_tg(context, mouse_2d)
        if tg_hit:
            return None  # Ray hits TubeGroom, don't hit surface
    ray_origin = view3d_utils.region_2d_to_origin_3d(context.region, context.region_data, mouse_2d)
    ray_direction = view3d_utils.region_2d_to_vector_3d(context.region, context.region_data, mouse_2d)
    if not ray_origin or not ray_direction:
        return None
    matrix_inv = obj.matrix_world.inverted()
    success, hit_loc, *_ = obj.ray_cast(
        matrix_inv @ ray_origin,
        (matrix_inv.to_3x3() @ ray_direction).normalized(),
    )
    return obj.matrix_world @ hit_loc if success else None
def ray_tg(context, mouse_2d):
    tg_object = geometry.get_tg_obj(context, allow_create=False)
    if not tg_object or not tg_object.data.vertices:
        return None, None
    ray_origin = view3d_utils.region_2d_to_origin_3d(context.region, context.region_data, mouse_2d)
    ray_direction = view3d_utils.region_2d_to_vector_3d(context.region, context.region_data, mouse_2d)
    if not ray_origin or not ray_direction:
        return None, None
    matrix_inv = tg_object.matrix_world.inverted()
    success, hit_loc, _, face_idx = tg_object.ray_cast(
        matrix_inv @ ray_origin,
        (matrix_inv.to_3x3() @ ray_direction).normalized()
    )
    if success:
        return tg_object, {'hit_loc': hit_loc, 'face_idx': face_idx}
    return None, None

# Selection utilities
def get_point(mouse_2d, context, distance=40, set_highlight=False):
    from . import operators
    distance_sq = distance * distance
    candidate = None
    min_dist_sq = float('inf')
    surface_obj = get_surface_obj(context)
    root_offset = 0.0
    
    # Get camera position for depth checking
    region = context.region
    region_data = context.region_data
    view_matrix = region_data.view_matrix
    camera_loc = view_matrix.inverted().translation
    
    # Get depth from raycast to TubeGroom mesh
    tg_object, hit_data = ray_tg(context, mouse_2d)
    hit_depth = None
    if tg_object and hit_data:
        hit_loc_world = tg_object.matrix_world @ hit_data['hit_loc']
        hit_depth = (hit_loc_world - camera_loc).length
    
    # Use 2D projection but respect depth - only select points at or in front of mesh depth
    for region_id, region in geometry.TubeGroom.regions.items():
        for subregion_id, subregion in region.subregions.items():
            if operators.modal_state.creation.current_region_points and subregion_id != 1:
                continue
            # When dragging, only allow snapping to root points to avoid self-intersection.
            if operators.modal_state.dragging_point and subregion_id != 1:
                continue
            for point_index, point in enumerate(subregion.points):
                if operators.modal_state.dragging_point and region_id == operators.modal_state.selection.region_id and subregion_id == operators.modal_state.selection.subregion_id and point_index == operators.modal_state.selection.point_index:
                    continue
                disp_pos = point.position
                if subregion_id == 1 and surface_obj and root_offset >= 0.0 and (hit_n := closest_point(surface_obj, point.position)) and hit_n[0] and hit_n[1]:
                    disp_pos = hit_n[0] + (hit_n[1] * root_offset)
                
                # Check depth - only consider points at or in front of mesh
                point_depth = (disp_pos - camera_loc).length
                if hit_depth is not None and point_depth > hit_depth + 0.02:  # 2cm tolerance behind mesh
                    continue
                    
                if (point_2d := project_2d(disp_pos, context)):
                    d_sq = (Vector(mouse_2d) - point_2d).length_squared
                    if d_sq < distance_sq and d_sq < min_dist_sq:
                        min_dist_sq = d_sq
                        candidate = {'rid': region_id, 'sid': subregion_id, 'pidx': point_index, 'pos3d': disp_pos}
    # Check current region points using 2D projection (since they don't have mesh yet)
    for j, point_3d in enumerate(operators.modal_state.creation.current_region_points):
        if operators.modal_state.dragging_point and operators.modal_state.selection.region_id == -1 and operators.modal_state.selection.point_index == j:
            continue
        disp_pos = point_3d
        if surface_obj and root_offset >= 0.0 and (hit_n := closest_point(surface_obj, point_3d)) and hit_n[0] and hit_n[1]:
            disp_pos = hit_n[0] + (hit_n[1] * root_offset)
        if (point_2d := project_2d(disp_pos, context)) and (d_sq := (Vector(mouse_2d) - point_2d).length_squared) < distance_sq and d_sq < min_dist_sq:
            min_dist_sq = d_sq
            candidate = {'rid': -1, 'sid': -1, 'pidx': j, 'pos3d': disp_pos}
    
    if not candidate:
        if set_highlight:
            operators.modal_state.creation.snap_to_nearest = None
        return -1, -1, -1, None

    if set_highlight:
        operators.modal_state.creation.snap_to_nearest = candidate['pos3d']

    return candidate['rid'], candidate['sid'], candidate['pidx'], candidate['pos3d']
def get_edge(mouse_2d, context, distance=50):
    tg_object, hit_data = ray_tg(context, mouse_2d)
    if not tg_object or not hit_data:
        return -1, -1, -1, None
    mesh = tg_object.data
    face = mesh.polygons[hit_data['face_idx']]
    hit_loc = hit_data['hit_loc']
    closest_edge_key = min(face.edge_keys, key=lambda ek: dist2_seg3(hit_loc, mesh.vertices[ek[0]].co, mesh.vertices[ek[1]].co), default=None)
    if not closest_edge_key:
        return -1, -1, -1, None
    # Project that 3D edge to 2D and check distance to the mouse.
    v1_world = tg_object.matrix_world @ mesh.vertices[closest_edge_key[0]].co
    v2_world = tg_object.matrix_world @ mesh.vertices[closest_edge_key[1]].co
    v1_2d = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, v1_world)
    v2_2d = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, v2_world)
    if not v1_2d or not v2_2d:
        return -1, -1, -1, None
    dist_2d_sq = dist2_seg2(mouse_2d, v1_2d, v2_2d)
    if dist_2d_sq > distance * distance:
        return -1, -1, -1, None
    # If close enough, extract the region/subregion/point IDs from the edge's vertices
    v_idx1, v_idx2 = closest_edge_key
    # Get attributes from both vertices
    rid1 = mesh.attributes['region_id'].data[v_idx1].value
    sid1 = mesh.attributes['subregion_id'].data[v_idx1].value
    rid2 = mesh.attributes['region_id'].data[v_idx2].value
    sid2 = mesh.attributes['subregion_id'].data[v_idx2].value
    if rid1 == rid2 and rid1 > 0 and sid1 != sid2 and sid1 > 0 and sid2 > 0:
        return rid1, min(sid1, sid2), max(sid1, sid2), 'vertical'
    if rid1 != rid2 or sid1 != sid2 or rid1 <= 0 or sid1 <= 0:
        return -1, -1, -1, None
    rid, sid = rid1, sid1
    region = geometry.TubeGroom.regions.get(rid)
    if not region:
        return -1, -1, -1, None
    subregion = region.subregions.get(sid)
    if not subregion:
        return -1, -1, -1, None
    v1_co = tg_object.matrix_world @ mesh.vertices[v_idx1].co
    v2_co = tg_object.matrix_world @ mesh.vertices[v_idx2].co
    pidx1, pidx2 = -1, -1
    for i, p in enumerate(subregion.points):
        if (p.position - v1_co).length_squared < 1e-8:
            pidx1 = i
        if (p.position - v2_co).length_squared < 1e-8:
            pidx2 = i
        if pidx1 != -1 and pidx2 != -1:
            break
    if pidx1 == -1 or pidx2 == -1:
        return -1, -1, -1, None
    n_points = len(subregion.points)
    edge_idx = min(pidx1, pidx2) if abs(pidx1 - pidx2) == 1 else (n_points - 1 if (pidx1 == 0 and pidx2 == n_points - 1) or (pidx2 == 0 and pidx1 == n_points - 1) else -1)
    return (rid, sid, edge_idx, 'horizontal') if edge_idx != -1 else (-1, -1, -1, None)
def get_region(mouse_2d, context):
    tg_object, hit_data = ray_tg(context, mouse_2d)
    if not tg_object or not hit_data:
        return -1, -1, -1
    mesh = tg_object.data
    face = mesh.polygons[hit_data['face_idx']]
    rid_attr_data = mesh.attributes['region_id'].data
    sid_attr_data = mesh.attributes['subregion_id'].data
    if face.vertices:
        rid = rid_attr_data[face.vertices[0]].value
        all_sids = {sid_attr_data[v_idx].value for v_idx in face.vertices}
        positive_sids = {s for s in all_sids if s > 0}
        if not positive_sids:
            # This is a purely interpolated face. We need to find the real subregions
            # by traversing the mesh edges.
            adj_list = [[] for _ in range(len(mesh.vertices))]
            for edge in mesh.edges:
                v1, v2 = edge.vertices
                adj_list[v1].append(v2)
                adj_list[v2].append(v1)
            found_sids = set()
            q = deque(face.vertices)
            visited = set(q)
            while q:
                curr_v_idx = q.popleft()
                curr_sid = sid_attr_data[curr_v_idx].value
                if curr_sid > 0:
                    found_sids.add(curr_sid)
                    if len(found_sids) >= 2:
                        break  # We found both boundaries
                    continue  # Don't traverse past a real subregion
                for neighbor_idx in adj_list[curr_v_idx]:
                    if neighbor_idx not in visited:
                        visited.add(neighbor_idx)
                        # Only add neighbors from the same region
                        if rid_attr_data[neighbor_idx].value == rid:
                            q.append(neighbor_idx)
            return (rid, min(found_sids), max(found_sids)) if len(found_sids) >= 2 else (-1, -1, -1)
        if len(positive_sids) == 1:
            sid = positive_sids.pop()
            return rid, sid, sid
        return rid, min(positive_sids), max(positive_sids)
    return -1, -1, -1
def get_region_collapsed(region_id, subregion_id, context, threshold=50):
    # Checks if a subregion's 2D bounding box is smaller than a pixel threshold.
    region = geometry.TubeGroom.regions.get(region_id)
    if not region:
        return False
    subregion = region.subregions.get(subregion_id)
    if not subregion:
        return False
    positions = subregion.get_positions()
    if not positions:
        return False
    # Project all points to 2D with single comprehension
    points_2d = [p2d for p in positions if (p2d := project_2d(p, context)) is not None]
    if len(points_2d) < 2:
        return True
    xs = [p.x for p in points_2d]
    ys = [p.y for p in points_2d]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    diagonal_sq = (max_x - min_x)**2 + (max_y - min_y)**2
    return diagonal_sq < threshold * threshold
def mouse_preview(context, event):
    """Manages visual previews for new points on the surface and snapping."""
    from . import operators
    from mathutils import Vector
    creation = operators.modal_state.creation
    
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)

    creation.temp_point = None
    creation.snap_target = None
    creation.snap_to_nearest = None
    
    current_pts = creation.current_region_points
    if len(current_pts) >= 3:
        first_2d = project_2d(current_pts[0], context)
        if first_2d and (Vector(mouse_2d) - Vector(first_2d)).length < 20:
            creation.snap_target = current_pts[0]
            creation.temp_point = Vector(current_pts[0])

    rid, sid, pidx, nearest_pos = get_point(mouse_2d, context, 50, True)
    if nearest_pos:
        creation.temp_point = Vector(nearest_pos)

    base_obj = get_surface_obj(context)
    if base_obj and base_obj.data.vertices:
        # Check if ray hits TubeGroom object first
        tg_hit, _ = ray_tg(context, mouse_2d)
        if not tg_hit:  # Only raycast to surface if no TubeGroom hit
            ray_origin = view3d_utils.region_2d_to_origin_3d(
                context.region, context.region_data, mouse_2d
            )
            ray_dir = view3d_utils.region_2d_to_vector_3d(
                context.region, context.region_data, mouse_2d
            )

            if ray_origin and ray_dir:
                inv = base_obj.matrix_world.inverted()
                success, hit_location, *_ = base_obj.ray_cast(
                    inv @ ray_origin,
                    (inv.to_3x3() @ ray_dir).normalized()
                )

                if success and creation.temp_point is None:
                    creation.temp_point = base_obj.matrix_world @ hit_location
                
    context.area.tag_redraw()

# Mesh attribute utilities
def update_attr(mesh, attr_name, attr_type, domain, values):
    """Create or update a mesh attribute."""
    if not values:
        return
    attrs = mesh.attributes
    if attr_name not in attrs:
        attrs.new(name=attr_name, type=attr_type, domain=domain)
    attr_data = attrs[attr_name].data
    if attr_type == 'FLOAT_COLOR':
        flat_values = [c for col in values for c in col]
        attr_data.foreach_set('color', flat_values)
    else:
        attr_data.foreach_set('value', values)
