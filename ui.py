import bpy
from mathutils import Vector
from bpy_extras import view3d_utils
from . import core as shared_data, geometry as geom, utils, interpolation, operators

# Operators
class TUBEGROOM_OT_create_tubegroom(bpy.types.Operator):
    bl_idname = "strand.create_tubegroom"
    bl_label = "Create TubeGroom"
    bl_description = "Create a new TubeGroom object for the selected mesh"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        target = context.scene.strand_raycast_target
        if not target or target.type != 'MESH' or utils.is_tubegroom_object(target):
            self.report({'ERROR'}, "Set target mesh")
            return {'CANCELLED'}
        if bpy.data.objects.get(f"GEO_TubeGroom_{target.name}"):
            self.report({'WARNING'}, "TubeGroom already exists")
            return {'CANCELLED'}
        obj = geom.get_or_create_tubegroom_object(context, target)
        if not obj:
            self.report({'ERROR'}, "Creation failed")
            return {'CANCELLED'}
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        context.view_layer.objects.active = obj
        self.report({'INFO'}, f"Created: {obj.name}")
        return {'FINISHED'}
class TUBEGROOM_OT_build_curves(bpy.types.Operator):
    bl_idname = "strand.build_curves"
    bl_label = "Build Curves"
    bl_description = "Generate the final hair curves object from the TubeGroom guides"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        if not context.scene.tubegroom_curves_enabled:
            self.report({'INFO'}, "Enable curves first")
            return {'CANCELLED'}
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        if not shared_data.regions:
            self.report({'WARNING'}, "No regions")
            return {'CANCELLED'}
        tg = context.active_object
        if not utils.is_tubegroom_object(tg):
            self.report({'ERROR'}, 'Select TubeGroom mesh')
            return {'CANCELLED'}
        base_name = utils.get_base_name_from_obj(tg)
        if not base_name:
            self.report({'ERROR'}, 'Invalid object name')
            return {'CANCELLED'}
        system = interpolation.generate_tubegroom_interpolation()
        if not system:
            self.report({'ERROR'}, 'Interpolation failed')
            return {'CANCELLED'}
        curves_obj = interpolation.build_curves_object_from_system(base_name, system)
        if not curves_obj:
            self.report({'ERROR'}, 'Build failed')
            return {'CANCELLED'}
        if hasattr(curves_obj.data, 'display') and hasattr(context.scene.render, 'hair_subdiv'):
            curves_obj.data.display.surface_subdivisions = context.scene.render.hair_subdiv
            if hasattr(curves_obj.data, 'update_tag'):
                curves_obj.data.update_tag()
        enable = context.scene.strand_interpolation_enabled
        curves_obj.show_in_front = enable
        if hasattr(curves_obj, 'color') and len(curves_obj.color) >= 4:
            curves_obj.color[3] = 0.8 if enable else 1.0
        self.report({'INFO'}, f'Built: {curves_obj.name}')
        return {'FINISHED'}

#-- Main Modal Operators --#
class TUBEGROOM_OT_edit_mode(bpy.types.Operator):
    bl_idname = "strand.edit_mode"
    bl_label = "Click to Add Points"
    bl_description = "Enter interactive edit mode to create and modify TubeGroom guides"
    bl_options = {'REGISTER', 'UNDO'}
    def modal(self, context, event):
        shared_data.just_cleared = False
        if shared_data.exit_request:
            shared_data.exit_request = False
            self.finish(context)
            context.area.tag_redraw()
            return {'FINISHED'}
        is_mouse_over_ui = any(r.type == 'UI' and r.x < event.mouse_x < r.x + r.width and r.y < event.mouse_y < r.y + r.height 
                               for r in context.area.regions)
        if self.mouse_was_over_ui and not is_mouse_over_ui:
            self.mouse_was_over_ui = False
            return {'PASS_THROUGH'}
        self.mouse_was_over_ui = is_mouse_over_ui
        if is_mouse_over_ui or context.region.type != 'WINDOW':
            if shared_data.temp_point:
                shared_data.temp_point = None
                context.area.tag_redraw()
            return {'PASS_THROUGH'}
        context.window.cursor_set('CROSSHAIR')
        
        if event.type in {'WHEELUPMOUSE', 'WHEELDOWNMOUSE', 'MIDDLEMOUSE'}:
            return {'PASS_THROUGH'}
        if event.type == 'MOUSEMOVE':
            if shared_data.extrude_mode:
                operators.handle_extrusion_mousemove(context, event, self.extrude_start_mouse_3d)
                context.area.tag_redraw()
            elif shared_data.scaling_mode:
                operators.handle_scaling_mousemove(context, event)
            elif shared_data.rotation_mode:
                operators.handle_rotation_mousemove(context, event)
            elif shared_data.dragging_point:
                operators.handle_vertex_drag(context, event)
            elif shared_data.move_subregion_mode:
                operators.handle_move_subregion_mousemove(context, event)
            elif shared_data.twist_mode:
                operators.handle_twist_mousemove(context, event)
            else:
                utils.handle_mouse_preview(context, event)
            context.area.tag_redraw()
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            mouse_2d = (event.mouse_region_x, event.mouse_region_y)
            shift, ctrl, alt = event.shift, event.ctrl, event.alt
            is_creating = bool(shared_data.current_region_points)
            rid_p, sid_p, pidx = utils.get_point_at_mouse_2d(
                mouse_2d,
                context,
                restrict_to_root_when_creating=is_creating,
            )
            rid_e, sid_e1, sid_e2_or_edge_idx, edge_type = utils.get_edge_at_mouse_2d_detailed(mouse_2d, context)
            rid_face, sid_face_min, sid_face_max = utils.get_region_at_mouse_2d(mouse_2d, context)
            body_sid = max(sid_face_max, sid_face_min) if rid_face > 0 and sid_face_min != sid_face_max else -1
            body_hit = body_sid > 1
            if ctrl and alt and not shift:
                if body_hit:
                    return self.handle_twist(context, event, rid_face, body_sid, hierarchical=True)
                if edge_type == 'horizontal':
                    if rid_e != -1 and sid_e1 > 1:
                        return self.handle_twist(context, event, rid_e, sid_e1)
                    if rid_e != -1 and sid_e1 == 1:
                        self.report({'WARNING'}, 'Cannot twist root')
                elif edge_type == 'vertical':
                    self.report({'WARNING'}, 'Use horizontal edges')
                else:
                    self.report({'WARNING'}, 'Cannot twist root')
                return {'RUNNING_MODAL'}
            if shift and alt and not ctrl:
                if pidx != -1 and sid_p > 0:
                    return self.start_move_column(context, event, rid_p, sid_p, pidx)
                if body_hit:
                    return self.handle_shift_alt_click(context, event, rid_face, body_sid, hierarchical=True)
                if edge_type == 'horizontal':
                    if sid_e1 > 1:
                        return self.handle_shift_alt_click(context, event, rid_e, sid_e1)
                    else:
                        self.report({'WARNING'}, 'Cannot rotate root')
                elif rid_e != -1:
                    self.report({'WARNING'}, 'Use horizontal edges')
                return {'RUNNING_MODAL'}

            if ctrl and not alt and not shift:
                return self.handle_ctrl_click(context, event)
            if alt and not ctrl and not shift:
                if body_hit:
                    return self.handle_alt_click(context, event, rid_face, body_sid, hierarchical=True)
                if pidx != -1 and sid_p > 1 and utils.is_subregion_small_on_screen(rid_p, sid_p, context, threshold=60):
                    return self.handle_alt_click(context, event, rid_p, sid_p)
                elif edge_type == 'horizontal':
                    if sid_e1 > 1:
                        return self.handle_alt_click(context, event, rid_e, sid_e1)
                    else:
                        self.report({'WARNING'}, 'Cannot scale root')
                elif rid_e != -1:
                    self.report({'WARNING'}, 'Use horizontal edges')
                return {'RUNNING_MODAL'}
            if shift and not ctrl and not alt:
                if body_hit:
                    return self.start_move_subregion(context, event, rid_face, body_sid, hierarchical=True)
                if pidx != -1:
                    is_small_point_subregion = (
                        sid_p > 1
                        and utils.is_subregion_small_on_screen(
                            rid_p,
                            sid_p,
                            context,
                            threshold=60,
                        )
                    )
                    if is_small_point_subregion:
                        return self.start_move_subregion(context, event, rid_p, sid_p)
                    else:
                        return self.start_drag(context, event, rid_p, sid_p, pidx)
                elif edge_type == 'horizontal' and sid_e1 != 1:
                    return self.start_move_subregion(context, event, rid_e, sid_e1)
                return {'RUNNING_MODAL'}
            if is_creating:
                if shared_data.snap_to_first and len(shared_data.current_region_points) >= 3:
                    return self.close_region(context)
                if shared_data.temp_point:
                    return self.add_point(context)
                return {'RUNNING_MODAL'}
            if rid_face > 0 and sid_face_min != sid_face_max:
                return self.handle_insert_subregion(context, rid_face, sid_face_min, sid_face_max)
            if rid_e != -1 and edge_type == 'horizontal':
                return self.handle_add_on_edge(context, rid_e, sid_e1, sid_e2_or_edge_idx)
            if shared_data.temp_point:
                return self.add_point(context)
            return {'RUNNING_MODAL'}
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            if shared_data.extrude_mode:
                return self.end_extrusion(context)
            elif shared_data.scaling_mode:
                return self.end_scaling(context)
            elif shared_data.rotation_mode:
                return self.end_rotation(context)
            elif shared_data.dragging_point:
                return self.end_drag(context)
            elif shared_data.move_subregion_mode:
                return self.end_move_subregion(context)
            elif shared_data.twist_mode:
                return self.end_twist(context)
        elif event.type == 'RIGHTMOUSE' and event.value == 'PRESS':
            result = operators.delete_item_at_mouse(context, event)
            if result:
                self.report({'INFO'}, result)
            return {'RUNNING_MODAL'}
        elif event.type == 'Z' and event.value == 'PRESS' and event.ctrl and not event.shift:
            return self.do_undo(context)
        elif event.type == 'Z' and event.value == 'PRESS' and event.ctrl and event.shift:
            return self.do_redo(context)
        elif event.type == 'ESC' and event.value == 'PRESS':
            return self.handle_escape(context)
        elif event.type == 'TAB' and event.value == 'PRESS':
            self.finish(context)
            return {'FINISHED'}
        # Pass through any unhandled events (system events, special inputs, etc.)
        return {'PASS_THROUGH'}

    #-- Modal Action Handlers --#
    def handle_add_on_edge(self, context, region_id, subregion_id, edge_index):
        ok = operators.add_hierarchical_point_on_edge(context, region_id, subregion_id, edge_index)
        self.report({'INFO' if ok else 'WARNING'}, 'Point inserted' if ok else 'Insert failed')
        return {'RUNNING_MODAL'}
    def handle_insert_subregion(self, context, region_id, subregion_id1, subregion_id2):
        ok = operators.insert_subregion_between(context, region_id, subregion_id1, subregion_id2, factor=0.5)
        self.report({'INFO' if ok else 'WARNING'}, 'Subregion inserted' if ok else 'Insert failed')
        return {'RUNNING_MODAL'}
    def handle_ctrl_click(self, context, event):
        ok, region_id, subregion_id, start_mouse_3d = operators.start_extrusion(context, event)
        if ok:
            self.extrude_start_mouse_3d = start_mouse_3d
            self.extrude_region_id = region_id
            self.report({'INFO'}, "Extrusion started")
        else:
            self.report({'WARNING'}, "No subregion to extrude")
        return {'RUNNING_MODAL'}
    def close_region(self, context):
        region_id = operators.end_region(context)
        if region_id is None:
            self.report({'WARNING'}, 'Need 3+ points')
            return {'RUNNING_MODAL'}
        self.report({'INFO'}, f"Region {region_id} closed")
        return {'RUNNING_MODAL'}
    def start_drag(self, context, event, region_id, subregion_id, point_index):
        if operators.start_drag(context, event, region_id, subregion_id, point_index):
            self.report({'INFO'}, "Dragging point")
        return {'RUNNING_MODAL'}
    def add_point(self, context):
        point_count = operators.add_point(context)
        self.report({'INFO'}, f"Point {point_count} added")
        return {'RUNNING_MODAL'}
    def end_drag(self, context):
        operators.end_drag(context)
        self.report({'INFO'}, "Moved")
        return {'RUNNING_MODAL'}
    def do_undo(self, context):
        return self._do_history_action(context, shared_data.restore_state, "Undone", "Nothing to undo")
    def do_redo(self, context):
        return self._do_history_action(context, shared_data.redo_state, "Redone", "Nothing to redo")
    def _do_history_action(self, context, action_func, success_msg, fail_msg):
        if action_func():
            geom.create_or_update_merged_mesh(context)
            from . import interpolation
            interpolation.update_tubegroom_interpolation(context, None, update_topology=True)
            context.area.tag_redraw()
            self.report({'INFO'}, success_msg)
        else:
            self.report({'INFO'}, fail_msg)
        return {'RUNNING_MODAL'}
    def handle_escape(self, context):
        result = operators.cancel_active_operation(context)
        if result:
            self.report({'INFO'}, result)
            context.area.tag_redraw()
        return {'RUNNING_MODAL'}
    def end_extrusion(self, context):
        operators.end_extrusion(context, self.extrude_region_id)
        self.report({'INFO'}, "Extrusion done")
        context.area.tag_redraw()
        return {'RUNNING_MODAL'}
    def start_move_column(self, context, event, region_id, subregion_id, point_index):
        if region_id == -1 or point_index == -1:
            return {'RUNNING_MODAL'}
        if operators.start_move_column(context, event, region_id, subregion_id, point_index):
            self.report({'INFO'}, 'Moving column')
        else:
            self.report({'WARNING'}, 'Cannot move column')
        return {'RUNNING_MODAL'}
    def start_move_subregion(self, context, event, region_id, subregion_id, hierarchical=False):
        region = shared_data.regions.get(region_id)
        if not region or subregion_id not in region.subregions:
            return {'RUNNING_MODAL'}
        positions = region.subregions[subregion_id].get_positions()
        center = sum(positions, Vector()) / len(positions) if positions else Vector()
        start_mouse_3d = view3d_utils.region_2d_to_location_3d(
            context.region, context.region_data,
            (event.mouse_region_x, event.mouse_region_y), center)
        if operators.start_move_subregion(region_id, subregion_id, start_mouse_3d, hierarchical=hierarchical):
            msg = 'Moving hierarchy' if hierarchical else 'Moving subregion'
            self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_move_subregion(self, context):
        operators.end_move_subregion(context)
        self.report({'INFO'}, 'Moved')
        return {'RUNNING_MODAL'}
    def handle_twist(self, context, event, region_id, subregion_id, hierarchical=False):
        shared_data.twist_mode = True
        shared_data.selected_region_id = region_id
        shared_data.selected_subregion_id = subregion_id
        shared_data.twist_start_mouse_x = event.mouse_region_x
        shared_data.twist_start_mouse_y = event.mouse_region_y
        operators.start_twist(context, region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Twisting hierarchy' if hierarchical else 'Twisting'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_twist(self, context):
        operators.end_twist(context)
        self.report({'INFO'}, 'Twist done')
        return {'RUNNING_MODAL'}
    def handle_alt_click(self, context, event, region_id, subregion_id, hierarchical=False):
        shared_data.scaling_mode = True
        shared_data.selected_region_id = region_id
        shared_data.selected_subregion_id = subregion_id
        shared_data.scaling_start_mouse_x = event.mouse_region_x
        shared_data.scaling_start_mouse_y = event.mouse_region_y
        self.scale_region_id = region_id
        self.scale_subregion_id = subregion_id
        operators.start_scaling(region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Scaling hierarchy' if hierarchical else 'Scaling'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_scaling(self, context):
        operators.end_scaling(context, self.scale_region_id)
        self.report({'INFO'}, "Scaling done")
        context.area.tag_redraw()
        return {'RUNNING_MODAL'}
    def handle_shift_alt_click(self, context, event, region_id, subregion_id, hierarchical=False):
        shared_data.rotation_mode = True
        shared_data.selected_region_id = region_id
        shared_data.selected_subregion_id = subregion_id
        shared_data.rotation_start_mouse_x = event.mouse_region_x
        self.rotation_region_id = region_id
        self.rotation_subregion_id = subregion_id
        operators.start_rotation(region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Rotating hierarchy' if hierarchical else 'Rotating'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}

    def end_rotation(self, context):
        operators.end_rotation(context, self.rotation_region_id, self.rotation_subregion_id)
        self.report({'INFO'}, "Rotation done")
        context.area.tag_redraw()
        return {'RUNNING_MODAL'}
    #-- Modal Invoke & Finish --#
    def invoke(self, context, event):
        if shared_data.edit_mode:
            shared_data.exit_request = True
            return {'FINISHED'}
        if context.area.type != 'VIEW_3D':
            return {'CANCELLED'}
        selected = context.active_object
        if not utils.is_tubegroom_object(selected):
            self.report({'ERROR'}, "Select TubeGroom object")
            return {'CANCELLED'}
        geom.ensure_region_attributes(selected.data)
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        from . import interpolation
        interpolation.rebuild_regions_from_merged_mesh(selected)
        surf_obj = utils.get_surface_object(context)
        if not surf_obj or surf_obj.type != 'MESH':
            self.report({'ERROR'}, "Select target object")
            return {'CANCELLED'}

        shared_data.root_offset = context.scene.strand_root_offset
        self._shading_backup = []
        areas = context.window.screen.areas if context.window else bpy.context.screen.areas
        for area in areas:
            if area.type != 'VIEW_3D':
                continue
            for space in area.spaces:
                if space.type != 'VIEW_3D' or not hasattr(space, 'shading'):
                    continue
                shade = space.shading
                self._shading_backup.append((space, shade.color_type, getattr(shade, 'color_attribute', None)))
                color_type_prop = shade.bl_rna.properties.get('color_type')
                allowed_types = {item.identifier for item in color_type_prop.enum_items} if color_type_prop else set()
                if 'VERTEX' in allowed_types:
                    shade.color_type = 'VERTEX'
                elif 'ATTRIBUTE' in allowed_types:
                    shade.color_type = 'ATTRIBUTE'
                    if hasattr(shade, 'color_attribute'):
                        shade.color_attribute = 'region_color'
        if hasattr(selected.data, 'color_attributes'):
            ca = selected.data.color_attributes
            if 'region_color' in ca:
                ca.active_color = ca['region_color']
        if context.mode == 'OBJECT':
            bpy.ops.object.select_all(action='DESELECT')
        # Initialize modal state variables.
        self.extrude_start_mouse_3d = None
        self.extrude_region_id = -1
        self.scale_region_id = -1
        self.scale_subregion_id = -1
        self.rotation_region_id = -1
        self.rotation_subregion_id = -1
        self.mouse_was_over_ui = False
        shared_data.edit_mode = True
        # Add drawing handlers for the UI.
        from . import drawing as draw_regions
        shared_data.regions_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_regions,
            (self, context),
            'WINDOW',
            'POST_VIEW',
        )
        shared_data.edges_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_edges,
            (self, context),
            'WINDOW',
            'POST_VIEW',
        )
        shared_data.points_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_points,
            (self, context),
            'WINDOW',
            'POST_VIEW',
        )
        context.window_manager.modal_handler_add(self)
        self.report({'INFO'}, "Edit mode")
        return {'RUNNING_MODAL'}
    def finish(self, context):
        for space, prev_type, prev_attr in self._shading_backup:
            if hasattr(space, 'shading'):
                shade = space.shading
                if prev_type is not None:
                    shade.color_type = prev_type
                if prev_attr is not None and hasattr(shade, 'color_attribute'):
                    shade.color_attribute = prev_attr
        self._shading_backup = []
        from . import drawing
        drawing.remove_draw_handlers()
        shared_data.clear_modal_state()
        shared_data.dragging_point = None
        from . import drawing as draw_regions
        draw_regions.clear_cache()
        shared_data.edit_mode = False
        if context.area:
            context.area.tag_redraw()
        total_points = len(shared_data.current_region_points)
        for region in shared_data.regions.values():
            for sub in region.subregions.values():
                total_points += len(sub.points)
        self.report({'INFO'}, f"{len(shared_data.regions)} regions, {total_points} points")
        return {'FINISHED'}
class TUBEGROOM_OT_mesh_to_tubegroom(bpy.types.Operator):
    bl_idname = "strand.convert"
    bl_label = "Mesh to TubeGroom"
    bl_description = "Convert a tube-like mesh into a TubeGroom region by detecting rings"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        from . import convert
        # Use the active object as source and the scene target as optional override
        src = context.active_object
        if not src or src.type != 'MESH':
            self.report({'ERROR'}, 'Select mesh to convert')
            return {'CANCELLED'}
        shared_data.reset_all_data(keep_edit_mode=False)
        target = getattr(context.scene, 'strand_raycast_target', None)
        ok, msg = convert.convert_mesh_object_to_tubegroom(context, src, target)
        if ok:
            if hasattr(src, 'hide_set'):
                src.hide_set(True)
            if hasattr(src, 'hide_viewport'):
                src.hide_viewport = True
            if hasattr(src, 'hide_render'):
                src.hide_render = True
            self.report({'INFO'}, msg)
            return {'FINISHED'}
        self.report({'ERROR'}, msg)
        return {'CANCELLED'}

class TUBEGROOM_OT_clear_points(bpy.types.Operator):
    bl_idname = "strand.clear_points"
    bl_label = "Clear All Regions"
    bl_description = "Delete all TubeGroom regions and clear all associated data"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        if not shared_data.edit_mode:
            from . import drawing
            drawing.remove_draw_handlers()
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        shared_data.save_state()
        shared_data.reset_all_data(keep_edit_mode=shared_data.edit_mode)

        for obj in bpy.data.objects:
            if not utils.is_tubegroom_object(obj):
                continue

            if obj.type == 'MESH' and obj.data:
                mesh = obj.data
                materials = list(mesh.materials)
                mesh.clear_geometry()
                for mat in materials:
                    if mat:
                        mesh.materials.append(mat)
                geom.ensure_region_attributes(mesh)

            base_name = utils.get_base_name_from_obj(obj)
            if base_name:
                cobj = bpy.data.objects.get(f"CRV_{base_name}")
                if cobj:
                    data_ref = cobj.data
                    bpy.data.objects.remove(cobj, do_unlink=True)
                    if data_ref and data_ref.users == 0 and hasattr(bpy.data, 'hair_curves'):
                        bpy.data.hair_curves.remove(data_ref)
        if context.area:
            context.area.tag_redraw()
        self.report({'INFO'}, "Cleared")
        return {'FINISHED'}

#-- UI Panels --#
class TUBEGROOM_PT_main_panel(bpy.types.Panel):
    bl_label = "TubeGroom"
    bl_idname = "TUBEGROOM_PT_main_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TubeGroom"

    
    def draw(self, context):
        layout = self.layout
        layout.label(text="TubeGroom Tools")
class TUBEGROOM_PT_mesh_panel(bpy.types.Panel):
    bl_label = "TubeGroom Mesh"
    bl_idname = "TUBEGROOM_PT_mesh_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TubeGroom"
    bl_parent_id = "TUBEGROOM_PT_main_panel"
    def draw(self, context):
        layout = self.layout

        box_create = layout.box()
        box_create.label(text="Target object")
        box_create.prop(context.scene, "strand_raycast_target", text="")
        box_create .operator("strand.create_tubegroom", text="Create TubeGroom", icon='OUTLINER_OB_MESH')
        box_create .operator('strand.convert', text='Mesh to TubeGroom', icon='MESH_CYLINDER')
        layout.separator()
        box_config = layout.box()
        box_config.label(text="Configuration:")
        box_config.prop(context.scene, "strand_root_offset_mm", text="Root Offset (mm)", slider=True)
        box_config.prop(context.scene, 'strand_view_segments', text="Resolution", slider=True)
        box_config.prop(context.scene, 'strand_bend_factor', text="Smooth", slider=True)
        box_ops = layout.box()
        # The edit button now acts as a toggle.
        col = box_ops.column(align=True)
        op_text = "Exit TubeGroom Edit" if shared_data.edit_mode else "TubeGroom Edit"
        op_icon = 'X' if shared_data.edit_mode else 'GREASEPENCIL'
        col .operator("strand.edit_mode", text=op_text, icon=op_icon)
        col .operator("strand.clear_points", text="Clear All Regions", icon='TRASH')

class TUBEGROOM_PT_interpolation_panel(bpy.types.Panel):
    bl_label = "TubeGroom Interpolation"
    bl_idname = "TUBEGROOM_PT_interpolation_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TubeGroom"
    bl_parent_id = "TUBEGROOM_PT_main_panel"
    bl_options = {'DEFAULT_CLOSED'}
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        box = layout.box()
        # Main toggle for the curves functionality
        box.prop(scene, 'tubegroom_curves_enabled', text='Enable Curves Display')
        # Sub-layout enabled when the main toggle is on
        sub_box = box.box()
        sub_box.active = scene.tubegroom_curves_enabled
        sub_box.prop(scene, 'strand_interpolation_enabled', text='Live Updates')
        sub_box.prop(scene, 'strand_density_factor', text='Density')
        # Shape settings collapsible section
        row = sub_box.row()
        row.alignment = 'LEFT' # Ensure the label is left-aligned
        row.prop(scene, "tubegroom_show_shape_settings", text=" Shape", # Added space for padding
                 icon='TRIA_RIGHT' if not scene.tubegroom_show_shape_settings else 'TRIA_DOWN',
                 emboss=False)
        if scene.tubegroom_show_shape_settings:
            shape_box = sub_box.box()
            # These controls mimic the old hair display settings for the curves object
            shape_box.label(text="Viewport Display")
            row = shape_box.row()
            if hasattr(scene.render, 'hair_type'):
                row.prop(scene.render, "hair_type", expand=True)
            if hasattr(scene.render, 'hair_subdiv'):
                shape_box.prop(scene.render, "hair_subdiv", text="Additional Subdivisions")
        sub_box .operator('strand.build_curves', text='Create Curves Object', icon='OUTLINER_DATA_CURVES')
class TUBEGROOM_PT_controls_panel(bpy.types.Panel):
    bl_label = "TubeGroom Controls"
    bl_idname = "TUBEGROOM_PT_controls_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TubeGroom"
    bl_parent_id = "TUBEGROOM_PT_main_panel"
    bl_options = {'DEFAULT_CLOSED'}
    def draw(self, _context):
        layout = self.layout
        box = layout.box()
        col = box.column(align=True)
        controls = [
            "Left Click: Add point",
            "Esc: Cancel Region",
            "Left Click: Insert",
            "Ctrl+Click: Extrude",
            "Shift+Click: Move (edge) / hierarchy (body)",
            "Alt+Click: Scale (edge) / hierarchy (body)",
            "Shift+Alt+Click: Rotate (edge/body) / Move column (point)",
            "Ctrl+Alt+Click: Twist (edge) / hierarchy (body)",
            "Right Click: Delete",
            "Ctrl+Z: Undo",
            "Ctrl+Shift+Z: Redo",
            "Tab: Exiting",
        ]
        for text in controls:
            command, _, description = text.partition(":")
            row = col.row(align=True)
            split = row.split(factor=0.45, align=True)
            split.label(text=command)
            split.label(text=description.strip())

# Registration
classes = [
    TUBEGROOM_PT_main_panel,
    TUBEGROOM_PT_mesh_panel,
    TUBEGROOM_PT_interpolation_panel,
    TUBEGROOM_PT_controls_panel,
    TUBEGROOM_OT_create_tubegroom,
    TUBEGROOM_OT_build_curves,
    TUBEGROOM_OT_edit_mode,
    TUBEGROOM_OT_mesh_to_tubegroom,
    TUBEGROOM_OT_clear_points,
]
def register():
    # Load custom icons
    import os
    import bpy.utils.previews
    if not shared_data.custom_icons:
        shared_data.custom_icons = bpy.utils.previews.new()
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    icon_path = os.path.join(icons_dir, "TG_Main_Icon.png")
    if os.path.exists(icon_path):
        shared_data.custom_icons.load("main_icon", icon_path, 'IMAGE')
    for cls in classes:
        bpy.utils.register_class(cls)
def unregister():
    shared_data.exit_request = True
    from . import drawing
    drawing.remove_draw_handlers()
    shared_data.clear_modal_state()

    # Remove custom icons preview collection
    import bpy.utils.previews
    if shared_data.custom_icons:
        bpy.utils.previews.remove(shared_data.custom_icons)
        shared_data.custom_icons = None
    for cls in reversed(classes):
        if hasattr(bpy.utils, 'unregister_class'):
            bpy.utils.unregister_class(cls)
