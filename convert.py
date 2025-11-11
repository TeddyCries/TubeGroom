import bmesh
from . import core as shared_data
from . import geometry as geom
from . import interpolation as _interp, utils

def propagate_rings_by_vertex_map_indices(bm, ordered_root_idx, allowed_idx=None):
    """Propagate rings using vertex indices within the same BMesh."""
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    allowed = None if allowed_idx is None else set(allowed_idx)
    current_ring_idx = [i for i in ordered_root_idx if (allowed is None or i in allowed)]
    n = len(current_ring_idx)
    if n < 3:
        return []

    rings_idx = [list(current_ring_idx)]

    def edge_between_idx(ia, ib):
        a = bm.verts[ia]
        b = bm.verts[ib]
        for e in a.link_edges:
            if b in e.verts:
                return e
        return None

    prev_faces = []
    for i in range(n):
        ia = current_ring_idx[i]
        ib = current_ring_idx[(i + 1) % n]
        e = edge_between_idx(ia, ib)
        if not e:
            return rings_idx
        faces = [f for f in e.link_faces
                 if len(f.verts) == 4 and (allowed is None or all(v.index in allowed for v in f.verts))]
        if not faces:
            return rings_idx
        prev_faces.append(faces[0])

    while True:
        next_ring_idx = []
        next_prev_faces = []
        boundary_hit = True
        for i in range(n):
            ia = current_ring_idx[i]
            ib = current_ring_idx[(i + 1) % n]
            e = edge_between_idx(ia, ib)
            if not e:
                return rings_idx
            faces = [f for f in e.link_faces
                     if len(f.verts) == 4 and (allowed is None or all(v.index in allowed for v in f.verts))]
            if not faces:
                return rings_idx
            fwd = faces[0] if len(faces) == 1 else (faces[1] if faces[0] is prev_faces[i] else faces[0])

            opp = utils.opposite_edge(fwd, e)
            if opp is None:
                return rings_idx
            if not opp.is_boundary:
                boundary_hit = False

            a = bm.verts[ia]
            b = bm.verts[ib]
            cand_v = None
            for ed in fwd.edges:
                if a in ed.verts:
                    w = ed.other_vert(a)
                    if w is not b:
                        cand_v = w
                        break
            if cand_v is None:
                return rings_idx

            next_ring_idx.append(cand_v.index)
            next_prev_faces.append(fwd)

        rings_idx.append(next_ring_idx)
        if boundary_hit:
            break
        current_ring_idx = next_ring_idx
        prev_faces = next_prev_faces

    return rings_idx

def choose_root_loop_by_distance(loops_idx_sets, bm, mesh_obj, target_obj):
    """Select the loop closest to the target surface using raycast."""
    bm.verts.ensure_lookup_table()
    
    if not target_obj or getattr(target_obj, 'type', None) != 'MESH':
        min_z = None
        root_set = None
        for idx_set in loops_idx_sets:
            avg_z = sum((mesh_obj.matrix_world @ bm.verts[i].co).z for i in idx_set) / len(idx_set)
            if min_z is None or avg_z < min_z:
                min_z = avg_z
                root_set = idx_set
        return root_set
    
    min_dist = None
    root_set = None
    target_matrix_inv = target_obj.matrix_world.inverted()
    
    for idx_set in loops_idx_sets:
        total_dist = 0
        for i in idx_set:
            co_world = mesh_obj.matrix_world @ bm.verts[i].co
            co_local = target_matrix_inv @ co_world
            success, closest_local, *_ = target_obj.closest_point_on_mesh(co_local)
            if success:
                closest_world = target_obj.matrix_world @ closest_local
                total_dist += (co_world - closest_world).length
        
        avg_dist = total_dist / len(idx_set)
        if min_dist is None or avg_dist < min_dist:
            min_dist = avg_dist
            root_set = idx_set
    
    return root_set

def convert_mesh_object_to_tubegroom(context, mesh_obj, target_obj=None):
    if not mesh_obj or getattr(mesh_obj, 'type', None) != 'MESH':
        return False, 'Select a mesh object to convert.'
    if len(mesh_obj.data.vertices) < 3:
        return False, 'Mesh has too few vertices to form rings.'

    islands = utils.find_islands(mesh_obj)
    all_region_rings = []

    bm = bmesh.new()
    try:
        bm.from_mesh(mesh_obj.data)
        bm.normal_update()
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()

        base_for_distance = getattr(context.scene, 'strand_raycast_target', None)

        for island in islands:
            loops_idx_sets = utils.find_boundary_loops(bm, allowed_idx=island)
            if not loops_idx_sets:
                continue
            root_idx_set = choose_root_loop_by_distance(loops_idx_sets, bm, mesh_obj, base_for_distance)
            if root_idx_set is None or len(root_idx_set) < 3:
                continue

            ordered_root_idx = utils.order_ring_from_indices(bm, root_idx_set)
            if not ordered_root_idx or len(ordered_root_idx) < 3:
                continue

            if base_for_distance and getattr(base_for_distance, 'type', None) == 'MESH':
                target_matrix_inv = base_for_distance.matrix_world.inverted()
                dmin = None
                start_i = 0
                for i, vi in enumerate(ordered_root_idx):
                    co_world = mesh_obj.matrix_world @ bm.verts[vi].co
                    co_local = target_matrix_inv @ co_world
                    success, closest_local, *_ = base_for_distance.closest_point_on_mesh(co_local)
                    if success:
                        closest_world = base_for_distance.matrix_world @ closest_local
                        d = (co_world - closest_world).length
                        if dmin is None or d < dmin:
                            dmin = d
                            start_i = i
                ordered_root_idx = ordered_root_idx[start_i:] + ordered_root_idx[:start_i]

            bm_rings_idx = propagate_rings_by_vertex_map_indices(bm, ordered_root_idx, allowed_idx=island)
            if not bm_rings_idx or len(bm_rings_idx) < 2:
                continue
            rings = [[mesh_obj.matrix_world @ bm.verts[i].co for i in ring] for ring in bm_rings_idx]
            all_region_rings.append(rings)
    finally:
        bm.free()

    if not all_region_rings:
        return False, 'Could not detect tube rings; ensure meshes are open tubes with quads.'
    
    tg_obj = geom.get_or_create_tubegroom_object(context, target_obj=base_for_distance, allow_create=True)
    if tg_obj is None:
        return False, 'Failed to create/find TubeGroom object.'

    created = 0
    for rings in all_region_rings:
        rid = shared_data.next_region_id
        region = shared_data.Region(rid, shared_data.current_region_color_index)
        for ring in rings:
            sub = region.create_subregion()
            for p in ring:
                sub.add_point(p.copy())
        shared_data.regions[rid] = region
        shared_data.next_region_id += 1
        shared_data.current_region_color_index += 1
        created += 1

    geom.create_or_update_merged_mesh(context)
    if getattr(context.scene, 'tubegroom_curves_enabled', True):
        base_name = utils.get_base_name_from_obj(tg_obj)
        system = _interp.generate_tubegroom_interpolation()
        if system and base_name:
            _interp.build_curves_object_from_system(base_name, system)
    return True, f'Created {created} regions from mesh.'
