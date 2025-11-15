import bmesh
from mathutils import Vector, Matrix
from bpy_extras import view3d_utils
from . import utils, geometry

# Global history for undo/redo
history = []
history_index = -1

class ModalState:
    def __init__(self):
        # Modal modes
        self.edit_mode = self.dragging_point = self.extrude_mode = False
        self.scaling_mode = self.rotation_mode = self.move_subregion_mode = False
        self.twist_mode = self.column_drag_mode = False

        # History NO está aquí - es global y persiste

        # Point selection state
        self.selected_region_id = self.selected_subregion_id = -1
        self.selected_point_index = -1
        self.current_region_points = []

        # Drag state
        self.drag_start_mouse_3d = None
        self.drag_original_position = None
        self.temp_point = None
        self.snap_target = None
        self.snap_to_first = False
        self.snap_to_nearest = None

        # Transform operations
        self.scaling_original_positions = None
        self.rotation_original_positions = None
        self.move_subregion_original_positions = None
        self.twist_original_positions = None
        self.twist_axis = None
        self.start_mouse_pos = None

        # Extrusion state
        self.extrude_preview_subregion_id = None
        self.extrude_normal = None
        self.extrude_height = 0.0

        # Column drag state
        self.column_drag_original_positions = None
        self.column_drag_active_subregion = -1

        # Color state
        self.current_region_color = None
        self.current_region_color_index = None

        # Mouse state
        self.mouse_was_over_ui = False
        
        # Highlight state
        self.highlight_region_id = -1
        self.highlight_subregion_id = -1

modal_state = ModalState()

# State management functions
def save_state():
    global history, history_index
    state = {
        'regions': {rid: {
            'region_id': r.region_id,
            'color_index': r.color_index,
            'next_subregion_id': r.next_subregion_id,
            'subregions': {sid: {
                'region_id': s.region_id,
                'subregion_id': s.subregion_id,
                'points': [p.position.copy() for p in s.points]
            } for sid, s in r.subregions.items()}
        } for rid, r in geometry.regions.items()},
        'next_region_id': geometry.next_region_id,
        'current_region_points': [p.copy() for p in modal_state.current_region_points]
    }
    history = history[:history_index + 1]
    history.append(state)
    history_index += 1
    if len(history) > 20:
        history.pop(0)
        history_index -= 1
def undo_state():
    global history, history_index
    if history_index <= 0 or not modal_state.edit_mode:
        return False
    history_index -= 1
    state = history[history_index]
    
    geometry.regions.clear()
    for rid, rdata in state['regions'].items():
        region = geometry.Region(rdata['region_id'], rdata['color_index'])
        region.next_subregion_id = rdata['next_subregion_id']
        for sid, sdata in rdata['subregions'].items():
            subregion = geometry.Subregion(sdata['region_id'], sdata['subregion_id'])
            for pos in sdata['points']:
                subregion.add_point(pos)
            region.subregions[sid] = subregion
        geometry.regions[rid] = region
    geometry.next_region_id = state['next_region_id']
    modal_state.current_region_points = [p.copy() for p in state.get('current_region_points', [])]
    clear_modal_state(keep_edit_mode=True)
   
    return True
def redo_state():
    global history, history_index
    if history_index >= len(history) - 1 or not modal_state.edit_mode:
        return False
    history_index += 1
    state = history[history_index]
    
    geometry.regions.clear()
    for rid, rdata in state['regions'].items():
        region = geometry.Region(rdata['region_id'], rdata['color_index'])
        region.next_subregion_id = rdata['next_subregion_id']
        for sid, sdata in rdata['subregions'].items():
            subregion = geometry.Subregion(sdata['region_id'], sdata['subregion_id'])
            for pos in sdata['points']:
                subregion.add_point(pos)
            region.subregions[sid] = subregion
        geometry.regions[rid] = region
    geometry.next_region_id = state['next_region_id']
    modal_state.current_region_points = [p.copy() for p in state.get('current_region_points', [])]
    clear_modal_state(keep_edit_mode=True)
    return True

# Clear modal state
def clear_modal_state(keep_edit_mode=None):
    global modal_state
    
    # Simply recreate the ModalState instance - clean and maintainable
    old_edit_mode = modal_state.edit_mode if keep_edit_mode is None else keep_edit_mode
    modal_state = ModalState()
    modal_state.edit_mode = old_edit_mode
def reset_all_data(keep_edit_mode=False, clear_history=False):
    global history, history_index
    geometry.regions.clear()
    geometry.next_region_id = 1
    # Simply use clear_modal_state - it recreates the instance cleanly
    clear_modal_state(keep_edit_mode=keep_edit_mode)
    if clear_history:
        history.clear()
        history_index = -1

# Region update functions
def update_geometry(context, region_id=None, update_topology=False):
    """Update mesh and interpolation for regions."""
    geometry.update_mesh_date(context, region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, region_id, update_topology=update_topology)
def update_all_geometry(context, update_topology=False):
    """Update merged mesh and interpolation for all regions."""
    update_geometry(context, update_topology=update_topology)

# Helper functions for region operations
def rebuild_points_with_insert(subregion, insert_index, position):
    ordered_points = [p.position.copy() for p in subregion.points]
    ordered_points.insert(insert_index + 1, position.copy())
    subregion.points.clear()
    for pos in ordered_points:
        subregion.add_point(pos)
    subregion.touch()
def collect_subregion_snapshots(region, start_sub_id, hierarchical, compute_pivots=False):
    snapshots = {}
    for sid in sorted(region.subregions.keys()):
        if sid < start_sub_id or (not hierarchical and sid != start_sub_id):
            continue
        sub = region.subregions[sid]
        points = [p.position.copy() for p in sub.points]
        if points:
            snapshots[sid] = points
    return snapshots, None

# Region creation operations
def add_point(context):
    if not modal_state.temp_point:
        return 0
    if not modal_state.current_region_points:
        modal_state.current_region_color_index = geometry.next_region_id
    modal_state.current_region_points.append(modal_state.temp_point.copy())
    modal_state.temp_point = None
    modal_state.snap_target = None
    context.area.tag_redraw()
    save_state()
    return len(modal_state.current_region_points)
def end_region(context):
    if len(modal_state.current_region_points) < 3:
        return None
    
    region_id = geometry.next_region_id
    region = geometry.Region(region_id, region_id)
    subregion = region.create_subregion()
    
    for point_pos in modal_state.current_region_points:
        subregion.add_point(point_pos)

    geometry.regions[region_id] = region
    geometry.next_region_id += 1
    
    # Usar clear_modal_state() para limpiar variables temporales
    clear_modal_state()
    modal_state.edit_mode = True
    geometry.update_mesh_date(context)
    context.area.tag_redraw()
    save_state()
    return region_id

# Insertion operations
def insert_point_edge(context, region_id, subregion_id, edge_index):
    """Insert a new point on an edge, propagating it to all subregions in the hierarchy."""
    region = geometry.regions.get(region_id)
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
    changed = False
    surface_obj = utils.get_surface_obj(context)
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
        geometry.update_mesh_date(context, region_id)
        from . import interpolation
        interpolation.update_tubegroom_interpolation(context, region_id, update_topology=True)
        save_state()
    return changed
def insert_subregion(context, region_id, subregion_id1, subregion_id2, factor=0.5):
    """Insert a new subregion between two existing ones."""
    region = geometry.regions.get(region_id)
    if not region:
        return False
    s1 = region.subregions.get(subregion_id1)
    s2 = region.subregions.get(subregion_id2)
    if not s1 or not s2 or len(s1.points) != len(s2.points) or len(s1.points) < 3:
        return False
    a = s1.get_positions()
    b = s2.get_positions()
    n_points = len(a)
    bend, _ = utils.get_bend_and_segs(context)
    num_to_insert = 1
    prev_sub = region.subregions.get(subregion_id1 - 1)
    next_sub = region.subregions.get(subregion_id2 + 1)
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
    new_sub = geometry.Subregion(region_id, new_sid)
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
    save_state()
    geometry.update_mesh_date(context, region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, region_id, update_topology=True)
    return True
def cancel_active_operation(context):
    if any([modal_state.dragging_point, modal_state.move_subregion_mode, modal_state.extrude_mode,
            modal_state.scaling_mode, modal_state.rotation_mode, modal_state.twist_mode]):
        undo_state()
        clear_modal_state(keep_edit_mode=True)
        for region_id in geometry.regions.keys():
            geometry.update_mesh_date(context, region_id)
        from . import interpolation
        interpolation.update_tubegroom_interpolation(context, None, update_topology=True)
        return "Operation cancelled"
    if modal_state.current_region_points:
        # Usar clear_modal_state() para limpiar variables temporales
        clear_modal_state()
        return "Region creation cancelled"
    return None

# Deletion operations
def delete_from_current_region(context, point_index):
    if len(modal_state.current_region_points) > 1:
        modal_state.current_region_points.pop(point_index)
    else:
        # Usar clear_modal_state() para limpiar variables temporales
        clear_modal_state()
        context.area.tag_redraw()
    save_state()
def delete_from_hierarchy(context, region_id, point_index):
    region = geometry.regions.get(region_id)
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
    if not region.subregions:
        del geometry.regions[region_id]
        msg = f"Region {region_id} deleted"
    else:
        msg = "Point deleted from hierarchy"
    update_all_geometry(context, update_topology=True)
    save_state()
    return msg
def delete_subregion(context, region_id, subregion_id):
    region = geometry.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return f"Subregion {region_id}-{subregion_id} not found."
    del region.subregions[subregion_id]
    region.touch()
    if not region.subregions:
        del geometry.regions[region_id]
    region.renumber_subregions()
    update_all_geometry(context, update_topology=True)
    context.area.tag_redraw()
    save_state()
    return f"Subregion {region_id}-{subregion_id} deleted."
def delete_item_at_mouse(context, event):
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    region_id_p, _, point_index_p, _ = utils.get_point(mouse_2d, context)
    if point_index_p != -1:
        result = delete_from_current_region(context, point_index_p) if region_id_p == -1 else delete_from_hierarchy(context, region_id_p, point_index_p)
        if result:
            context.area.tag_redraw()
        return result
    region_id_e, subregion_id_e, _, edge_type = utils.get_edge(mouse_2d, context)
    if edge_type == 'horizontal' and region_id_e > 0:
        return delete_subregion(context, region_id_e, subregion_id_e) if subregion_id_e > 1 else delete_region(context, region_id_e)
    return None
def delete_region(context, region_id):
    if region_id in geometry.regions:
        del geometry.regions[region_id]
        update_all_geometry(context, update_topology=True)
        save_state()
        return f"Region {region_id} deleted"
    return None

# Vertex drag operation
def start_drag(context, event, region_id, subregion_id, point_index):
    if region_id == -1 or point_index == -1:
        return False
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.selected_point_index = point_index
    modal_state.dragging_point = True
    original_pos = geometry.regions[region_id].subregions[subregion_id].points[point_index].position
    initial_3d = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), original_pos)
    if not initial_3d:
        return False
    modal_state.drag_start_mouse_3d = initial_3d
    modal_state.drag_original_position = original_pos.copy()
    context.area.tag_redraw()
    return True
def start_move_column(context, event, region_id, subregion_id, point_index):
    region = geometry.regions.get(region_id)
    if not region or subregion_id not in region.subregions or point_index < 0:
        return False
    column_orig = {sid: sub.points[point_index].position.copy() 
                   for sid in sorted(region.subregions.keys()) 
                   if sid >= subregion_id and point_index < len((sub := region.subregions[sid]).points)}
    if not column_orig:
        return False
    base_pos = column_orig.get(subregion_id, next(iter(column_orig.values())))
    modal_state.column_drag_mode = True
    modal_state.column_drag_original_positions = column_orig
    modal_state.column_drag_active_subregion = subregion_id
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.selected_point_index = point_index
    modal_state.dragging_point = True
    modal_state.drag_original_position = base_pos.copy()
    modal_state.drag_start_mouse_3d = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), base_pos)
    context.area.tag_redraw()
    return True
def handle_vertex_drag(context, event):
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    if modal_state.column_drag_mode:
        if not modal_state.column_drag_original_positions or modal_state.selected_region_id not in geometry.regions or modal_state.selected_point_index < 0:
            return
        region = geometry.regions[modal_state.selected_region_id]
        base_original = modal_state.column_drag_original_positions.get(modal_state.column_drag_active_subregion, modal_state.drag_original_position)
        new_position = modal_state.snap_to_nearest or (
            utils.ray_surface(context, mouse_2d, utils.get_surface_obj(context), avoid_tg=False) if modal_state.column_drag_active_subregion == 1 else
            base_original + (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, modal_state.drag_original_position or Vector((0, 0, 0))) - modal_state.drag_start_mouse_3d)
            if modal_state.drag_start_mouse_3d else None
        )
        if new_position:
            displacement = new_position - base_original
            for sid, orig in modal_state.column_drag_original_positions.items():
                if sid in region.subregions and modal_state.selected_point_index < len(region.subregions[sid].points):
                    region.subregions[sid].points[modal_state.selected_point_index].position = orig + displacement
                    region.subregions[sid].touch()
            update_geometry(context, modal_state.selected_region_id)
            context.area.tag_redraw()
        return
    utils.get_point(mouse_2d, context, 15, True)
    subregion = geometry.regions[modal_state.selected_region_id].subregions[modal_state.selected_subregion_id]
    if 0 <= modal_state.selected_point_index < len(subregion.points):
        new_position = modal_state.snap_to_nearest or (
            utils.ray_surface(context, mouse_2d, utils.get_surface_obj(context), avoid_tg=False) if modal_state.selected_subregion_id == 1 else
            modal_state.drag_original_position + (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, modal_state.drag_original_position or Vector((0, 0, 0))) - modal_state.drag_start_mouse_3d)
            if modal_state.drag_start_mouse_3d else None
        )
        if new_position:
            subregion.points[modal_state.selected_point_index].position = new_position
            subregion.touch()
            update_geometry(context, modal_state.selected_region_id)
            context.area.tag_redraw()
def end_drag(context):
    region_id = modal_state.selected_region_id
    if region_id in geometry.regions:
        update_geometry(context, region_id)
    save_state()
    clear_modal_state()
    context.area.tag_redraw()
    return True

# Move subregion operation
def start_move_subregion(region_id, subregion_id, start_mouse_3d, hierarchical=False):
    region = geometry.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    modal_state.move_subregion_mode = True
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.move_subregion_original_positions = snapshots
    modal_state.drag_start_mouse_3d = start_mouse_3d
    modal_state.drag_original_position = sum(snapshots[subregion_id], Vector()) / len(snapshots[subregion_id])
    return True
def handle_move_subregion_mousemove(context, event):
    if not modal_state.move_subregion_mode or modal_state.selected_region_id not in geometry.regions or not modal_state.move_subregion_original_positions:
        return False
    displacement = view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), modal_state.drag_original_position or Vector((0, 0, 0))) - modal_state.drag_start_mouse_3d
    region = geometry.regions[modal_state.selected_region_id]
    for sid, originals in modal_state.move_subregion_original_positions.items():
        if sid in region.subregions:
            for i, orig_pos in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = orig_pos + displacement
            region.subregions[sid].touch()
    geometry.update_mesh_date(context, modal_state.selected_region_id)
    from . import interpolation
    interpolation.update_tubegroom_interpolation(context, modal_state.selected_region_id)
    return True
def end_move_subregion(context):
    if modal_state.selected_region_id in geometry.regions:
        update_geometry(context, modal_state.selected_region_id)
    save_state()
    clear_modal_state()
    return True

# Extrusion operation
def start_extrusion(context, event):
    mouse_2d = (event.mouse_region_x, event.mouse_region_y)
    region_id, sid_min, sid_max = utils.get_region(mouse_2d, context)
    if region_id <= 0 or sid_min != sid_max or region_id not in geometry.regions:
        return False, None, None, None
    region_obj = geometry.regions[region_id]
    if not region_obj.subregions or sid_min != max(region_obj.subregions.keys()):
        return False, None, None, None
    tip_subregion = region_obj.subregions[sid_min]
    base_positions = tip_subregion.get_positions()
    if not base_positions:
        return False, None, None, None
    modal_state.extrude_normal = utils.tip_normal(base_positions, utils.get_surface_obj(context))
    preview_subregion = region_obj.create_subregion()
    for p in tip_subregion.points:
        preview_subregion.add_point(p.position.copy())
    modal_state.extrude_preview_subregion_id = preview_subregion.subregion_id
    modal_state.extrude_mode = True
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = sid_min
    modal_state.extrude_height = 0.0
    return True, region_id, sid_min, view3d_utils.region_2d_to_location_3d(context.region, context.region_data, mouse_2d, sum(base_positions, Vector()) / len(base_positions))
def handle_extrusion_mousemove(context, event, start_mouse_3d):
    if not modal_state.extrude_mode or not start_mouse_3d or not modal_state.extrude_normal or modal_state.extrude_normal.length == 0:
        return False
    modal_state.extrude_height = max(0.0, (view3d_utils.region_2d_to_location_3d(context.region, context.region_data, (event.mouse_region_x, event.mouse_region_y), start_mouse_3d) - start_mouse_3d).dot(modal_state.extrude_normal))
    if modal_state.selected_region_id in geometry.regions and modal_state.extrude_preview_subregion_id:
        region = geometry.regions[modal_state.selected_region_id]
        if modal_state.selected_subregion_id in region.subregions and modal_state.extrude_preview_subregion_id in region.subregions:
            tip_subregion = region.subregions[modal_state.selected_subregion_id]
            preview_sub = region.subregions[modal_state.extrude_preview_subregion_id]
            for i, tip_point in enumerate(tip_subregion.points):
                if i < len(preview_sub.points):
                    preview_sub.points[i].position = tip_point.position + (modal_state.extrude_normal * modal_state.extrude_height)
            preview_sub.touch()
            update_geometry(context, modal_state.selected_region_id, update_topology=True)
            return True
    return False
def end_extrusion(context, region_id):
    if region_id in geometry.regions:
        update_geometry(context, region_id, update_topology=True)
    save_state()
    clear_modal_state()
    return True

# Scaling operation
def start_scaling(region_id, subregion_id, hierarchical=False):
    region = geometry.regions.get(region_id)
    if not region or subregion_id not in region.subregions or subregion_id == 1:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    modal_state.scaling_mode = True
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.scaling_original_positions = snapshots
    # start_mouse_pos will be set in handle_scaling_mousemove
    return True
def handle_scaling_mousemove(context, event):
    if not modal_state.scaling_mode or modal_state.selected_region_id not in geometry.regions or not modal_state.scaling_original_positions:
        return False
    if modal_state.start_mouse_pos is None:
        modal_state.start_mouse_pos = (event.mouse_region_x, event.mouse_region_y)
        return True
    
    dx = event.mouse_region_x - modal_state.start_mouse_pos[0]
    dy = event.mouse_region_y - modal_state.start_mouse_pos[1]
    scale_factor = max(0.001, 1.0 + (dx + dy) * 0.005)
    
    region = geometry.regions[modal_state.selected_region_id]
    for sid, originals in modal_state.scaling_original_positions.items():
        if sid in region.subregions:
            pivot = sum(originals, Vector()) / len(originals)  # Calculate pivot dynamically
            for i, initial_pos in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = pivot + (initial_pos - pivot) * scale_factor
            region.subregions[sid].touch()
    update_geometry(context, modal_state.selected_region_id)
    return True
def end_scaling(context, region_id):
    if region_id in geometry.regions:
        update_geometry(context, region_id)
    save_state()
    clear_modal_state()
    return True

# Rotation operation
def start_rotation(region_id, subregion_id, hierarchical=True):
    region = geometry.regions.get(region_id)
    if not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    modal_state.rotation_mode = True
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.rotation_original_positions = snapshots
    # start_mouse_pos will be set in handle_rotation_mousemove
    return True
def handle_rotation_mousemove(context, event):
    if not modal_state.rotation_mode or modal_state.selected_region_id not in geometry.regions or not modal_state.rotation_original_positions:
        return False
    if modal_state.start_mouse_pos is None:
        modal_state.start_mouse_pos = (event.mouse_region_x, event.mouse_region_y)
        return True
    
    angle = (event.mouse_region_x - modal_state.start_mouse_pos[0]) * -0.01
    pivot = sum(modal_state.rotation_original_positions[modal_state.selected_subregion_id], Vector()) / len(modal_state.rotation_original_positions[modal_state.selected_subregion_id])  # Calculate pivot dynamically
    rot_m = Matrix.Identity(4)
    if abs(angle) > 0.001:
        view_axis = context.region_data.view_matrix.inverted().to_3x3() @ Vector((0, 0, 1))
        rot_m = Matrix.Rotation(angle, 4, view_axis)
    
    region = geometry.regions[modal_state.selected_region_id]
    for sub_id, originals in modal_state.rotation_original_positions.items():
        if sub_id in region.subregions:
            for i, orig in enumerate(originals):
                if i < len(region.subregions[sub_id].points):
                    region.subregions[sub_id].points[i].position = pivot + (rot_m @ (orig - pivot)).to_3d()
            region.subregions[sub_id].touch()
    update_geometry(context, modal_state.selected_region_id)
    return True
def end_rotation(context, region_id, subregion_id):
    if region_id in geometry.regions:
        update_geometry(context, region_id)
    save_state()
    clear_modal_state()
    return True

# Twist operation
def start_twist(context, region_id, subregion_id, hierarchical=True):
    region = geometry.regions.get(region_id)
    if subregion_id == 1 or not region or subregion_id not in region.subregions:
        return False
    snapshots, _ = collect_subregion_snapshots(region, subregion_id, hierarchical)
    if subregion_id not in snapshots:
        return False
    modal_state.twist_mode = True
    modal_state.selected_region_id = region_id
    modal_state.selected_subregion_id = subregion_id
    modal_state.twist_axis = utils.tip_normal(snapshots[subregion_id], utils.get_surface_obj(context))
    modal_state.twist_original_positions = snapshots
    # start_mouse_pos will be set in handle_twist_mousemove
    return True
def handle_twist_mousemove(context, event):
    if not modal_state.twist_mode or modal_state.selected_region_id not in geometry.regions or not modal_state.twist_original_positions:
        return False
    if modal_state.start_mouse_pos is None:
        modal_state.start_mouse_pos = (event.mouse_region_x, event.mouse_region_y)
        return True
    
    dx = event.mouse_region_x - modal_state.start_mouse_pos[0]
    dy = event.mouse_region_y - modal_state.start_mouse_pos[1]
    angle = (dx + dy) * -0.002
    
    anchor = sum(modal_state.twist_original_positions[modal_state.selected_subregion_id], Vector()) / len(modal_state.twist_original_positions[modal_state.selected_subregion_id])  # Calculate anchor dynamically
    rot_m = Matrix.Identity(4)
    if abs(angle) > 0.001 and modal_state.twist_axis and modal_state.twist_axis.length > 0:
        rot_m = Matrix.Rotation(angle, 4, modal_state.twist_axis.normalized())
    
    region = geometry.regions[modal_state.selected_region_id]
    for sid, originals in modal_state.twist_original_positions.items():
        if sid in region.subregions:
            for i, orig in enumerate(originals):
                if i < len(region.subregions[sid].points):
                    region.subregions[sid].points[i].position = anchor + (rot_m @ (orig - anchor)).to_3d()
            region.subregions[sid].touch()
    update_geometry(context, modal_state.selected_region_id)
    return True
def end_twist(context):
    if modal_state.selected_region_id in geometry.regions:
        update_geometry(context, modal_state.selected_region_id)
    save_state()
    clear_modal_state()
    return True

# Mesh conversion operation
def propagate_rings_by_vertex_map_indices(bm, ordered_root_idx, allowed_idx=None):
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
        a, b = bm.verts[ia], bm.verts[ib]
        for e in a.link_edges:
            if b in e.verts:
                return e
        return None

    prev_faces = []
    for i in range(n):
        ia, ib = current_ring_idx[i], current_ring_idx[(i + 1) % n]
        e = edge_between_idx(ia, ib)
        if not e:
            return rings_idx
        faces = [f for f in e.link_faces if len(f.verts) == 4 and (allowed is None or all(v.index in allowed for v in f.verts))]
        if not faces:
            return rings_idx
        prev_faces.append(faces[0])

    while True:
        next_ring_idx = []
        next_prev_faces = []
        boundary_hit = True
        for i in range(n):
            ia, ib = current_ring_idx[i], current_ring_idx[(i + 1) % n]
            e = edge_between_idx(ia, ib)
            if not e:
                return rings_idx
            faces = [f for f in e.link_faces if len(f.verts) == 4 and (allowed is None or all(v.index in allowed for v in f.verts))]
            if not faces:
                return rings_idx
            fwd = faces[0] if len(faces) == 1 else (faces[1] if faces[0] is prev_faces[i] else faces[0])

            opp = utils.opposite_edge(fwd, e)
            if opp is None or opp.is_boundary:
                boundary_hit = boundary_hit and opp is not None
            else:
                boundary_hit = False

            a, b = bm.verts[ia], bm.verts[ib]
            cand_v = next((ed.other_vert(a) for ed in fwd.edges if a in ed.verts and (cand_v := ed.other_vert(a)) is not b), None)
            if not cand_v:
                return rings_idx

            next_ring_idx.append(cand_v.index)
            next_prev_faces.append(fwd)

        rings_idx.append(next_ring_idx)
        if boundary_hit:
            break
        current_ring_idx, prev_faces = next_ring_idx, next_prev_faces

    return rings_idx
def choose_root_loop_by_distance(loops_idx_sets, bm, mesh_obj, target_obj):
    bm.verts.ensure_lookup_table()
    min_dist = root_set = None

    for idx_set in loops_idx_sets:
        total_dist = sum((mesh_obj.matrix_world @ bm.verts[i].co - target_obj.closest_point_on_mesh(mesh_obj.matrix_world @ bm.verts[i].co)[1]).length
                        for i in idx_set)
        avg_dist = total_dist / len(idx_set)
        if min_dist is None or avg_dist < min_dist:
            min_dist, root_set = avg_dist, idx_set

    return root_set
def convert_mesh_object_to_tubegroom(context, mesh_obj, target_obj=None):
    if not mesh_obj or mesh_obj.type != 'MESH' or len(mesh_obj.data.vertices) < 3:
        return False, 'Select a mesh object with at least 3 vertices to convert.'

    islands = utils.find_islands(mesh_obj)
    all_region_rings = []

    bm = bmesh.new()
    try:
        bm.from_mesh(mesh_obj.data)
        bm.normal_update()
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        base_for_distance = context.scene.strand_raycast_target

        for island in islands:
            loops_idx_sets = utils.boundary_loops(bm, allowed_idx=island)
            if not loops_idx_sets:
                continue
            root_idx_set = choose_root_loop_by_distance(loops_idx_sets, bm, mesh_obj, base_for_distance)
            if not root_idx_set or len(root_idx_set) < 3:
                continue

            ordered_root_idx = utils.order_ring_from_indices(bm, root_idx_set)
            if not ordered_root_idx or len(ordered_root_idx) < 3:
                continue

            if base_for_distance and base_for_distance.type == 'MESH':
                target_matrix_inv = base_for_distance.matrix_world.inverted()
                dmin = start_i = None
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
                if start_i:
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
    
    tg_obj = geometry.get_tg_obj(context, target_obj=base_for_distance, allow_create=True)
    if not tg_obj:
        return False, 'Failed to create/find TubeGroom object.'

    for rings in all_region_rings:
        rid = geometry.next_region_id
        region = geometry.Region(rid, rid)
        for ring in rings:
            sub = region.create_subregion()
            for p in ring:
                sub.add_point(p.copy())
        geometry.regions[rid] = region
        geometry.next_region_id += 1

    geometry.update_mesh_date(context)
    if context.scene.tubegroom_curves_enabled:
        base_name = utils.get_base_name(tg_obj)
        from . import interpolation
        system = interpolation.generate_tubegroom_interpolation()
        if system and base_name:
            interpolation.build_curves_object_from_system(base_name, system)
    return True, f'Created {len(all_region_rings)} regions from mesh.'