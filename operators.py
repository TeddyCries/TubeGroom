from mathutils import Vector
from bpy_extras import view3d_utils
from . import core as shared_data, utils, geometry

# CRUD Operations
def _update_all_geometry(context, update_topology=False):
    """Update merged mesh and interpolation for all regions."""
    geometry.create_or_update_merged_mesh(context)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, None, update_topology=update_topology)

def create_region_from_current():
    """Finalize current points into a new region."""
    if len(shared_data.current_region_points) < 3:
        return None
    region = shared_data.Region(shared_data.next_region_id, shared_data.current_region_color_index)
    subregion = region.create_subregion()
    for position in shared_data.current_region_points:
        subregion.add_point(position)
    region.touch()
    shared_data.regions[shared_data.next_region_id] = region
    result_id = shared_data.next_region_id
    shared_data.next_region_id += 1
    shared_data.current_region_color_index += 1
    shared_data.current_region_points.clear()
    return result_id

def add_point(context):
    """Add the current preview point to the active region being created."""
    save_state()
    if shared_data.temp_point is None:
        return 0
    shared_data.current_region_points.append(shared_data.temp_point)
    context.area.tag_redraw()
    return len(shared_data.current_region_points)

def end_region(context):
    """Finalize the current list of points into a new, permanent guide region."""
    if len(shared_data.current_region_points) >= 3:
        save_state()
        region_id = create_region_from_current()
        shared_data.temp_point = None
        shared_data.snap_to_first = False
        if region_id is not None:
            shared_data.edit_mode = True
            geometry.create_or_update_merged_mesh(context)
        context.area.tag_redraw()
        return region_id
    return None

def add_hierarchical_point_on_edge(context, region_id, subregion_id, edge_index):
    """Insert a new point on an edge, propagating it to all subregions in the hierarchy."""
    region = shared_data.regions.get(region_id)
    if not region:
        return False
    base = region.subregions.get(subregion_id)
    if not base:
        return False
    base_positions = base.get_positions()
    n = len(base_positions)
    if n < 3:
        return False
    i = int(edge_index) % n
    j = (i + 1) % n
    save_state()
    changed = False
    surface_obj = utils.get_surface_object(context)
    for sid in sorted(region.subregions.keys()):
        sub = region.subregions[sid]
        pos = sub.get_positions()
        if len(pos) != n:
            continue
        p_mid = (pos[i] + pos[j]) * 0.5
        p_new = p_mid
        if sid == 1 and surface_obj and surface_obj.type == 'MESH':
            matrix_inv = surface_obj.matrix_world.inverted()
            local_point = matrix_inv @ p_mid
            success, hit_location_local, _, _ = surface_obj.closest_point_on_mesh(local_point)
            if success:
                p_new = surface_obj.matrix_world @ hit_location_local
        rebuild_points_with_insert(sub, i, p_new)
        changed = True
        sub.touch()
    if changed:
        geometry.update_region_mesh(context, region_id)
        from . import interpolation
        interpolation.update_tubegroom_interpolation(context, region_id, update_topology=True)
    return changed

def insert_subregion_between(context, region_id, subregion_id1, subregion_id2, factor=0.5):
    """Insert a new subregion between two existing ones."""
    region = shared_data.regions.get(region_id)
    if not region:
        return False
    s1 = region.subregions.get(subregion_id1)
    s2 = region.subregions.get(subregion_id2)
    if not s1 or not s2 or len(s1.points) != len(s2.points) or len(s1.points) < 3:
        return False
    a = s1.get_positions()
    b = s2.get_positions()
    n_points = len(a)
    bend, _ = utils.get_view_settings(context)
    num_to_insert = 1
    prev_sub = region.subregions.get(subregion_id1 - 1)
    next_sub = region.subregions.get(subregion_id2 + 1)
    save_state()
    insert_sid = max(subregion_id1, subregion_id2)
    new_subregions = {}
    for sid, sub in region.subregions.items():
        if sid < insert_sid:
            new_subregions[sid] = sub
    subs_to_move = sorted([sub for sid, sub in region.subregions.items() if sid >= insert_sid], key=lambda s: s.subregion_id)
    for sub in subs_to_move:
        new_id = sub.subregion_id + num_to_insert
        sub.subregion_id = new_id
        for p in sub.points:
            p.subregion_id = new_id
        new_subregions[new_id] = sub
    region.subregions = new_subregions
    p_m = prev_sub.get_positions() if prev_sub and len(prev_sub.points) == n_points else a
    p_n = next_sub.get_positions() if next_sub and len(next_sub.points) == n_points else b
    new_sid = insert_sid
    new_sub = shared_data.Subregion(region_id, new_sid)
    t = factor
    if bend == 0.0:
        for j in range(n_points):
            new_sub.add_point(a[j].lerp(b[j], t))
    else:
        t2, t3 = t * t, t * t * t
        h00, h10 = 2.0 * t3 - 3.0 * t2 + 1.0, t3 - 2.0 * t2 + t
        h01, h11 = -2.0 * t3 + 3.0 * t2, t3 - t2
        for j in range(n_points):
            p_a, p_b = a[j], b[j]
            p_a_prev = p_m[j] if j < len(p_m) else p_a
            p_b_next = p_n[j] if j < len(p_n) else p_b
            t_a, t_b = (p_b - p_a_prev) * 0.5, (p_b_next - p_a) * 0.5
            p_hermite = (p_a * h00) + (t_a * h10) + (p_b * h01) + (t_b * h11)
            p_lin = p_a.lerp(p_b, t)
            new_sub.add_point(p_lin.lerp(p_hermite, bend))
    new_sub.touch()
    region.subregions[new_sid] = new_sub
    region.next_subregion_id += num_to_insert
    geometry.update_region_mesh(context, region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, region_id, update_topology=True)
    return True

def rebuild_points_with_insert(subregion, insert_index, position):
    """Helper to rebuild a subregion's points list after inserting a new point."""
    ordered_points = [p.position.copy() for p in subregion.points]
    ordered_points.insert(insert_index + 1, position.copy())
    subregion.points.clear()
    for pos in ordered_points:
        subregion.add_point(pos)
    subregion.touch()

def delete_from_current_region(point_index):
    """Delete a point from the guide region currently being created."""
    if len(shared_data.current_region_points) > 1:
        shared_data.current_region_points.pop(point_index)
    else:
        shared_data.current_region_points.clear()
        shared_data.temp_point = None
        shared_data.snap_to_first = False

def delete_from_hierarchy(context, region_id, point_index):
    """Delete a point and its corresponding points throughout the entire region hierarchy."""
    region = shared_data.regions.get(region_id)
    if not region:
        return
    for sid in list(region.subregions.keys()):
        if 0 <= point_index < len(region.subregions[sid].points):
            del region.subregions[sid].points[point_index]
            region.subregions[sid].touch()
    sids_to_delete = [sid for sid, sub in region.subregions.items() if len(sub.points) < 3]
    if sids_to_delete:
        for sid in sids_to_delete:
            del region.subregions[sid]
        region.touch()
        region.renumber_subregions()
        msg = "Point and subregions deleted"
    else:
        msg = "Point deleted from hierarchy"
    if not region.subregions:
        del shared_data.regions[region_id]
        msg = f"Region {region_id} deleted"
    _update_all_geometry(context, update_topology=True)
    return msg

def delete_subregion(context, region_id, subregion_id):
    """Delete a specific subregion by its ID."""
    region = shared_data.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return f"Subregion {region_id}-{subregion_id} not found."
    del region.subregions[subregion_id]
    region.touch()
    if not region.subregions:
        del shared_data.regions[region_id]
    region.renumber_subregions()
    _update_all_geometry(context, update_topology=True)
    context.area.tag_redraw()
    return f"Subregion {region_id}-{subregion_id} deleted."

def delete_item_at_mouse(context, event):
    """Main deletion logic: decides what to delete based on what's under the mouse."""
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    region_id_p, _, point_index_p = utils.get_point_at_mouse_2d(mouse_2d, context)
    if point_index_p != -1:
        save_state()
        if region_id_p == -1:
            result = delete_from_current_region(point_index_p)
        else:
            result = delete_from_hierarchy(context, region_id_p, point_index_p)
        if result:
            context.area.tag_redraw()
        return result
    region_id_e, subregion_id_e, _, edge_type = utils.get_edge_at_mouse_2d_detailed(mouse_2d, context)
    if edge_type == 'horizontal' and region_id_e > 0:
        save_state()
        return delete_subregion(context, region_id_e, subregion_id_e) if subregion_id_e > 1 else delete_region(context, region_id_e)
    return None

def delete_region(context, region_id):
    """Delete an entire guide region by its ID."""
    if region_id in shared_data.regions:
        del shared_data.regions[region_id]
        _update_all_geometry(context, update_topology=True)
        return f"Region {region_id} deleted"
    return None

# Helper functions
def _raycast_to_surface(context, mouse_2d, obj):
    """Perform raycast from mouse to surface object."""
    if not obj:
        return None
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

def _update_region_geometry(context, region_id, update_topology=False):
    """Update mesh and interpolation for a region."""
    geometry.update_region_mesh(context, region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, region_id, update_topology=update_topology)

# State & Operators Logic
def save_state():
    shared_data.save_state()

def cancel_active_operation(context):
    if (shared_data.dragging_point or shared_data.move_subregion_mode or
        shared_data.extrude_mode or shared_data.scaling_mode or
        shared_data.rotation_mode or shared_data.twist_mode):
        shared_data.restore_state()
        shared_data.clear_modal_state(keep_edit_mode=True)
        for region_id in shared_data.regions.keys():
            geometry.update_region_mesh(context, region_id)
        from . import interpolation
        interpolation.update_tubegroom_interpolation(context, None, update_topology=True)
        return "Operation cancelled"
    if shared_data.current_region_points:
        shared_data.current_region_points.clear()
        shared_data.temp_point = None
        shared_data.snap_to_first = False
        return "Region creation cancelled"
    return None

# Tool Handlers: Point Dragging
def start_drag(context, event, region_id, subregion_id, point_index):
    if region_id == -1 or point_index == -1:
        return False
    save_state()
    shared_data.selected_region_id = region_id
    shared_data.selected_subregion_id = subregion_id
    shared_data.selected_point_index = point_index
    shared_data.dragging_point = True
    original_pos = shared_data.regions[region_id].subregions[subregion_id].points[point_index].position
    initial_3d = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), original_pos)
    if not initial_3d:
        return False
    shared_data.drag_start_mouse_3d = initial_3d
    shared_data.drag_original_position = original_pos.copy()
    context.area.tag_redraw()
    return True

def start_move_column(context, event, region_id, subregion_id, point_index):
    region = shared_data.regions.get(region_id)
    if not region or subregion_id not in region.subregions or point_index < 0:
        return False
    column_orig = {sid: sub.points[point_index].position.copy() 
                   for sid in sorted(region.subregions.keys()) 
                   if sid >= subregion_id and point_index < len((sub := region.subregions[sid]).points)}
    if not column_orig:
        return False
    save_state()
    base_pos = column_orig.get(subregion_id, next(iter(column_orig.values())))
    shared_data.column_drag_mode = True
    shared_data.column_drag_original_positions = column_orig
    shared_data.column_drag_active_subregion = subregion_id
    shared_data.selected_region_id = region_id
    shared_data.selected_subregion_id = subregion_id
    shared_data.selected_point_index = point_index
    shared_data.dragging_point = True
    shared_data.drag_original_position = base_pos.copy()
    shared_data.drag_start_mouse_3d = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), base_pos)
    context.area.tag_redraw()
    return True

def handle_vertex_drag(context, event):
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    if getattr(shared_data, 'column_drag_mode', False):
        if not shared_data.column_drag_original_positions or shared_data.selected_region_id not in shared_data.regions or shared_data.selected_point_index < 0:
            return
        region = shared_data.regions[shared_data.selected_region_id]
        base_original = shared_data.column_drag_original_positions.get(shared_data.column_drag_active_subregion, shared_data.drag_original_position)
        new_position = shared_data.snap_to_existing or (
            _raycast_to_surface(context, mouse_2d, utils.get_surface_object(context)) if shared_data.column_drag_active_subregion == 1 else
            base_original + (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, shared_data.drag_original_position or Vector((0, 0, 0))) - shared_data.drag_start_mouse_3d)
            if shared_data.drag_start_mouse_3d else None
        )
        if new_position:
            displacement = new_position - base_original
            for sid, orig in shared_data.column_drag_original_positions.items():
                if sid in region.subregions and shared_data.selected_point_index < len(region.subregions[sid].points):
                    region.subregions[sid].points[shared_data.selected_point_index].position = orig + displacement
                    region.subregions[sid].touch()
            _update_region_geometry(context, shared_data.selected_region_id)
            context.area.tag_redraw()
        return
    utils.get_point_at_mouse_2d(mouse_2d, context, 15, True)
    subregion = shared_data.regions[shared_data.selected_region_id].subregions[shared_data.selected_subregion_id]
    if 0 <= shared_data.selected_point_index < len(subregion.points):
        new_position = shared_data.snap_to_existing or (
            _raycast_to_surface(context, mouse_2d, utils.get_surface_object(context)) if shared_data.selected_subregion_id == 1 else
            shared_data.drag_original_position + (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, shared_data.drag_original_position or Vector((0, 0, 0))) - shared_data.drag_start_mouse_3d)
            if shared_data.drag_start_mouse_3d else None
        )
        if new_position:
            subregion.points[shared_data.selected_point_index].position = new_position
            subregion.touch()
            _update_region_geometry(context, shared_data.selected_region_id)
            context.area.tag_redraw()

def end_drag(context):
    region_id = shared_data.selected_region_id
    if region_id in shared_data.regions:
        _update_region_geometry(context, region_id)
    shared_data.dragging_point = False
    shared_data.column_drag_mode = False
    shared_data.selected_region_id = -1
    shared_data.selected_subregion_id = -1
    shared_data.selected_point_index = -1
    shared_data.drag_start_mouse_3d = None
    shared_data.drag_original_position = None
    shared_data.column_drag_original_positions = None
    shared_data.column_drag_active_subregion = -1
    context.area.tag_redraw()
    return True

# Tool handlers - subregion move
def collect_subregion_snapshots(region, start_sub_id, hierarchical, compute_pivots=False):
    """Collect original positions (and optional pivots) for the target subregions."""
    snapshots = {}
    pivots = {} if compute_pivots else None
    for sid in sorted(region.subregions.keys()):
        if sid < start_sub_id:
            continue
        if not hierarchical and sid != start_sub_id:
            continue
        sub = region.subregions[sid]
        points = [p.position.copy() for p in sub.points]
        if not points:
            continue
        snapshots[sid] = points
        if compute_pivots:
            pivots[sid] = sum(points, Vector()) / len(points)
    return snapshots, pivots

def start_move_subregion(region_id, subregion_id, start_mouse_3d, hierarchical=False):
    region = shared_data.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    save_state()
    shared_data.move_subregion_mode = True
    shared_data.move_subregion_region_id = region_id
    shared_data.selected_subregion_id = subregion_id
    shared_data.selected_region_id = region_id
    shared_data.move_subregion_original_positions = snapshots
    shared_data.drag_start_mouse_3d = start_mouse_3d
    shared_data.drag_original_position = sum(snapshots[subregion_id], Vector()) / len(snapshots[subregion_id])
    return True

def handle_move_subregion_mousemove(context, event):
    if not shared_data.move_subregion_mode or shared_data.move_subregion_region_id not in shared_data.regions or not shared_data.move_subregion_original_positions:
        return False
    displacement = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), shared_data.drag_original_position or Vector((0, 0, 0))) - shared_data.drag_start_mouse_3d
    region = shared_data.regions[shared_data.move_subregion_region_id]
    for sid, originals in shared_data.move_subregion_original_positions.items():
        if sid in region.subregions:
            for i, orig_pos in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = orig_pos + displacement
            region.subregions[sid].touch()
    geometry.update_region_mesh(context, shared_data.move_subregion_region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, shared_data.move_subregion_region_id)
    return True

def end_move_subregion(context):
    if shared_data.move_subregion_region_id in shared_data.regions:
        _update_region_geometry(context, shared_data.move_subregion_region_id)
    shared_data.move_subregion_mode = False
    shared_data.move_subregion_region_id = -1
    shared_data.move_subregion_original_positions = None
    shared_data.drag_start_mouse_3d = None
    shared_data.drag_original_position = None
    return True

# Tool handlers - extrusion
def start_extrusion(context, event):
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    region_id, sid_min, sid_max = utils.get_region_at_mouse_2d(mouse_2d, context)
    if region_id <= 0 or sid_min != sid_max or region_id not in shared_data.regions:
        return False, None, None, None
    region_obj = shared_data.regions[region_id]
    if not region_obj.subregions or sid_min != max(region_obj.subregions.keys()):
        return False, None, None, None
    tip_subregion = region_obj.subregions[sid_min]
    base_positions = tip_subregion.get_positions()
    if not base_positions:
        return False, None, None, None
    save_state()
    shared_data.extrude_normal = utils.calculate_tip_face_normal(base_positions, utils.get_surface_object(context))
    preview_subregion = region_obj.create_subregion()
    for p in tip_subregion.points:
        preview_subregion.add_point(p.position.copy())
    shared_data.extrude_preview_subregion_id = preview_subregion.subregion_id
    shared_data.extrude_mode = True
    shared_data.selected_region_id = region_id
    shared_data.selected_subregion_id = sid_min
    shared_data.temp_extrude_height = 0.0
    return True, region_id, sid_min, view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, sum(base_positions, Vector()) / len(base_positions))

def handle_extrusion_mousemove(context, event, start_mouse_3d):
    if not shared_data.extrude_mode or not start_mouse_3d or not shared_data.extrude_normal or shared_data.extrude_normal.length == 0:
        return False
    extrude_height = max(0.0, (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), start_mouse_3d) - start_mouse_3d).dot(shared_data.extrude_normal))
    shared_data.temp_extrude_height = extrude_height
    if shared_data.selected_region_id in shared_data.regions and shared_data.extrude_preview_subregion_id:
        region = shared_data.regions[shared_data.selected_region_id]
        if shared_data.selected_subregion_id in region.subregions and shared_data.extrude_preview_subregion_id in region.subregions:
            tip_subregion = region.subregions[shared_data.selected_subregion_id]
            preview_sub = region.subregions[shared_data.extrude_preview_subregion_id]
            for i, tip_point in enumerate(tip_subregion.points):
                if i < len(preview_sub.points):
                    preview_sub.points[i].position = tip_point.position + (shared_data.extrude_normal * extrude_height)
            preview_sub.touch()
            _update_region_geometry(context, shared_data.selected_region_id, update_topology=True)
            return True
    return False

def end_extrusion(context, region_id):
    if region_id in shared_data.regions:
        _update_region_geometry(context, region_id, update_topology=True)
    shared_data.clear_modal_state()
    return True

# Tool handlers - scaling
def start_scaling(region_id, subregion_id, hierarchical=False):
    region = shared_data.regions.get(region_id)
    if not region or subregion_id not in region.subregions or subregion_id == 1:
        return False
    snapshots, pivots = collect_subregion_snapshots(region, subregion_id, hierarchical, compute_pivots=True)
    if subregion_id not in snapshots:
        return False
    save_state()
    shared_data.scaling_pivot_centers = pivots or {}
    shared_data.scaling_original_positions = snapshots
    return True

def handle_scaling_mousemove(context, event):
    if not shared_data.scaling_mode or shared_data.selected_region_id not in shared_data.regions:
        return False
    if not hasattr(shared_data, 'scaling_start_mouse_x'):
        shared_data.scaling_start_mouse_x = float(event.mouse_region_x)
        shared_data.scaling_start_mouse_y = float(event.mouse_region_y)
    scale_factor = max(0.001, 1.0 + ((float(event.mouse_region_x) - shared_data.scaling_start_mouse_x + float(event.mouse_region_y) - shared_data.scaling_start_mouse_y) * 0.005))
    region = shared_data.regions[shared_data.selected_region_id]
    for sid, originals in shared_data.scaling_original_positions.items():
        if sid in region.subregions and (pivot := shared_data.scaling_pivot_centers.get(sid)):
            for i, initial_pos in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = pivot + (initial_pos - pivot) * scale_factor
            region.subregions[sid].touch()
    _update_region_geometry(context, shared_data.selected_region_id)
    return True

def end_scaling(context, region_id):
    if region_id in shared_data.regions:
        _update_region_geometry(context, region_id)
    shared_data.clear_modal_state()
    return True

# Tool handlers - rotation and twist
def start_rotation(region_id, subregion_id, hierarchical=True):
    region = shared_data.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    save_state()
    shared_data.rotation_pivot_center = sum(snapshots[subregion_id], Vector()) / len(snapshots[subregion_id])
    shared_data.rotation_original_positions = snapshots
    return True

def handle_rotation_mousemove(context, event):
    if not shared_data.rotation_mode or shared_data.selected_region_id not in shared_data.regions:
        return False
    if not hasattr(shared_data, 'rotation_start_mouse_x'):
        shared_data.rotation_start_mouse_x = float(event.mouse_region_x)
    angle = (float(event.mouse_region_x) - shared_data.rotation_start_mouse_x) * -0.01
    from mathutils import Matrix
    rot_m = Matrix.Rotation(angle, 4, context.region_data.view_matrix.inverted().to_3x3() @ Vector((0, 0, 1))) if abs(angle) > 0.001 else Matrix.Identity(4)
    region = shared_data.regions[shared_data.selected_region_id]
    for sub_id, originals in shared_data.rotation_original_positions.items():
        if sub_id in region.subregions:
            for i, orig in enumerate(originals):
                if i < len(region.subregions[sub_id].points):
                    region.subregions[sub_id].points[i].position = shared_data.rotation_pivot_center + (rot_m @ (orig - shared_data.rotation_pivot_center)).to_3d()
            region.subregions[sub_id].touch()
    _update_region_geometry(context, shared_data.selected_region_id)
    return True

def end_rotation(context, region_id, subregion_id):
    if region_id in shared_data.regions:
        _update_region_geometry(context, region_id)
    shared_data.clear_modal_state()
    return True

def start_twist(context, region_id, subregion_id, hierarchical=True):
    region = shared_data.regions.get(region_id)
    if subregion_id == 1 or not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    save_state()
    shared_data.twist_mode = True
    shared_data.selected_region_id = region_id
    shared_data.selected_subregion_id = subregion_id
    shared_data.twist_anchor = sum(snapshots[subregion_id], Vector()) / len(snapshots[subregion_id])
    shared_data.twist_axis = utils.calculate_tip_face_normal(snapshots[subregion_id], utils.get_surface_object(context))
    shared_data.twist_original_positions = snapshots
    return True

def handle_twist_mousemove(context, event):
    if not shared_data.twist_mode or shared_data.selected_region_id not in shared_data.regions or shared_data.selected_subregion_id not in shared_data.regions[shared_data.selected_region_id].subregions:
        return False
    if not hasattr(shared_data, 'twist_start_mouse_x'):
        shared_data.twist_start_mouse_x = float(event.mouse_region_x)
        shared_data.twist_start_mouse_y = float(event.mouse_region_y)
    angle = ((float(event.mouse_region_x) - shared_data.twist_start_mouse_x + float(event.mouse_region_y) - shared_data.twist_start_mouse_y) * -0.005)
    from mathutils import Matrix
    rot_m = Matrix.Rotation(angle, 4, shared_data.twist_axis.normalized()) if abs(angle) > 0.0005 and shared_data.twist_axis and shared_data.twist_axis.length > 0 else Matrix.Identity(4)
    region = shared_data.regions[shared_data.selected_region_id]
    for sid, originals in shared_data.twist_original_positions.items():
        if sid in region.subregions:
            for i, orig in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = shared_data.twist_anchor + (rot_m @ (orig - shared_data.twist_anchor)).to_3d()
            region.subregions[sid].touch()
    _update_region_geometry(context, shared_data.selected_region_id)
    return True

def end_twist(context):
    if shared_data.selected_region_id in shared_data.regions:
        _update_region_geometry(context, shared_data.selected_region_id)
    shared_data.clear_modal_state()
    return True