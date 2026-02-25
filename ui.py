import bpy
from mathutils import Vector
from bpy_extras import view3d_utils
from . import geometry as geom, utils, operators

class TUBEGROOM_PT_main_panel(bpy.types.Panel):
    bl_label = "TubeGroom"
    bl_idname = "TUBEGROOM_PT_main_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "TubeGroom"

    def draw_header(self, context):
        self.layout.label(text="", icon_value=bpy.utils.previews.custom_icons["main_icon"].icon_id)
    
    def draw(self, context):
        pass   
class TUBEGROOM_OT_create_tubegroom(bpy.types.Operator):
    bl_idname = "strand.create_tubegroom"
    bl_label = "Create TubeGroom"
    bl_description = "Create a new TubeGroom object for the selected mesh"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        target = context.scene.strand_raycast_target
        if not target or target.type != 'MESH' or utils.tubegroom_object(target):
            self.report({'ERROR'}, "Set target mesh")
            return {'CANCELLED'}
        if bpy.data.objects.get(f"GEO_TubeGroom_{target.name}"):
            self.report({'WARNING'}, "TubeGroom already exists")
            return {'CANCELLED'}
        obj = geom.get_tg_obj(context, target)
        if not obj:
            self.report({'ERROR'}, "Creation failed")
            return {'CANCELLED'}
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        context.view_layer.objects.active = obj
        self.report({'INFO'}, f"Created: {obj.name}")
        return {'FINISHED'}
class TUBEGROOM_OT_mesh_to_tubegroom(bpy.types.Operator):
    bl_idname = "strand.convert"
    bl_label = "Mesh to TubeGroom"
    bl_description = "Convert a tube-like mesh into a TubeGroom region by detecting rings"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        src = context.active_object
        if not src or src.type != 'MESH':
            self.report({'ERROR'}, 'Select mesh to convert')
            return {'CANCELLED'}
        operators.reset_all_data(keep_edit_mode=False)
        target = context.scene.strand_raycast_target
        ok, msg = operators.convert_mesh_object_to_tubegroom(context, src, target)
        if ok:
            src.hide_set(True)
            src.hide_viewport = True
            src.hide_render = True
            self.report({'INFO'}, msg)
            return {'FINISHED'}
        self.report({'ERROR'}, msg)
        return {'CANCELLED'}

# Edit Mode Operator
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
        edit_mode = operators.modal_state.edit_mode
        op_text = "Exit TubeGroom Edit" if edit_mode else "TubeGroom Edit"
        op_icon = 'X' if edit_mode else 'GREASEPENCIL'
        col.operator("strand.edit_mode", text=op_text, icon=op_icon)
        col.operator("strand.clear_points", text="Clear All Regions", icon='TRASH')
class TUBEGROOM_OT_edit_mode(bpy.types.Operator):
    bl_idname = "strand.edit_mode"
    bl_label = "Click to Add Points"
    bl_description = "Enter interactive edit mode to create and modify TubeGroom guides"
    bl_options = {'REGISTER', 'UNDO'}

    def invoke(self, context, event):
        return self.start_modal(context, event)
    
    def modal(self, context, event):

        # Check if edit mode was disabled externally
        if not operators.modal_state.edit_mode:
            self.end_modal(context)
            context.area.tag_redraw()
            return {'FINISHED'}
        is_mouse_over_ui = any(r.type == 'UI' and r.x < event.mouse_x < r.x + r.width and r.y < event.mouse_y < r.y + r.height 
                               for r in context.area.regions)
        if operators.modal_state.mouse_was_over_ui and not is_mouse_over_ui:
            operators.modal_state.mouse_was_over_ui = False
            return {'PASS_THROUGH'}
        operators.modal_state.mouse_was_over_ui = is_mouse_over_ui
        if is_mouse_over_ui or context.region.type != 'WINDOW':
            if operators.modal_state.creation.temp_point:
                context.area.tag_redraw()
            return {'PASS_THROUGH'}
        context.window.cursor_set('CROSSHAIR')
        
        if event.type in {'WHEELUPMOUSE', 'WHEELDOWNMOUSE', 'MIDDLEMOUSE'}:
            return {'PASS_THROUGH'}
        if event.type == 'MOUSEMOVE':
            if operators.modal_state.extrude_mode:
                operators.handle_extrusion_mousemove(context, event, operators.modal_state.drag.start_mouse_3d)
                context.area.tag_redraw()
            elif operators.modal_state.scaling_mode:
                operators.handle_scaling_mousemove(context, event)
            elif operators.modal_state.rotation_mode:
                operators.handle_rotation_mousemove(context, event)
            elif operators.modal_state.dragging_point:
                operators.handle_vertex_drag(context, event)
                utils.mouse_preview(context, event)
            elif operators.modal_state.move_subregion_mode:
                operators.handle_move_subregion_mousemove(context, event)
            elif operators.modal_state.twist_mode:
                operators.handle_twist_mousemove(context, event)
            else:
                utils.mouse_preview(context, event)
            context.area.tag_redraw()
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            mouse_2d = (event.mouse_region_x, event.mouse_region_y)
            shift, ctrl, alt = event.shift, event.ctrl, event.alt
            is_creating = bool(operators.modal_state.creation.current_region_points)
            rid_p, sid_p, pidx, _ = utils.get_point(
                mouse_2d,
                context,
            )
            rid_e, sid_e1, sid_e2_or_edge_idx, edge_type = utils.get_edge(mouse_2d, context)
            rid_face, sid_face_min, sid_face_max = utils.get_region(mouse_2d, context)
            body_sid = max(sid_face_max, sid_face_min) if rid_face > 0 and sid_face_min != sid_face_max else -1
            body_hit = body_sid > 1
            if ctrl and alt and not shift:
                if body_hit:
                    return self.start_twist(context, event, rid_face, body_sid, hierarchical=True)
                if edge_type == 'horizontal':
                    if rid_e != -1 and sid_e1 > 1:
                        return self.start_twist(context, event, rid_e, sid_e1)
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
                    return self.start_rotation(context, event, rid_face, body_sid, hierarchical=True)
                if edge_type == 'horizontal':
                    if sid_e1 > 1:
                        return self.start_rotation(context, event, rid_e, sid_e1)
                    else:
                        self.report({'WARNING'}, 'Cannot rotate root')
                elif rid_e != -1:
                    self.report({'WARNING'}, 'Use horizontal edges')
                return {'RUNNING_MODAL'}

            if ctrl and not alt and not shift:
                return self.start_extrusion(context, event)
            if alt and not ctrl and not shift:
                if body_hit:
                    return self.start_scaling(context, event, rid_face, body_sid, hierarchical=True)
                if pidx != -1 and sid_p > 1 and utils.get_region_collapsed(rid_p, sid_p, context, threshold=60):
                    return self.start_scaling(context, event, rid_p, sid_p)
                elif edge_type == 'horizontal':
                    if sid_e1 > 1:
                        return self.start_scaling(context, event, rid_e, sid_e1)
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
                        and utils.get_region_collapsed(
                            rid_p,
                            sid_p,
                            context,
                            threshold=60,
                        )
                    )
                    if is_small_point_subregion:
                        return self.start_move_subregion(context, event, rid_p, sid_p)
                    else:
                        return self.start_drag_point(context, event, rid_p, sid_p, pidx)
                elif edge_type == 'horizontal' and sid_e1 != 1:
                    return self.start_move_subregion(context, event, rid_e, sid_e1)
                return {'RUNNING_MODAL'}
            
            # Double-click or simple click on edge to insert point
            if edge_type == 'horizontal' and rid_e != -1 and sid_e1 > 0:
                if operators.insert_point_edge(context, rid_e, sid_e1, sid_e2_or_edge_idx):
                    self.report({'INFO'}, 'Point inserted')
                    context.area.tag_redraw()
                    return {'RUNNING_MODAL'}
            
            # Double-click or simple click on body to insert subregion
            if body_hit and rid_face != -1 and rid_face in geom.TubeGroom.regions:
                adjacent_subregions = [sid for sid in geom.TubeGroom.regions[rid_face].subregions.keys() if sid != body_sid]
                if len(adjacent_subregions) >= 1:
                    closest_sub = min(adjacent_subregions, key=lambda s: abs(s - body_sid))
                    if operators.insert_subregion(context, rid_face, min(body_sid, closest_sub), max(body_sid, closest_sub)):
                        self.report({'INFO'}, 'Subregion inserted')
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
            
            if is_creating:
                if operators.modal_state.creation.snap_target and len(operators.modal_state.creation.current_region_points) >= 3:
                    return self.end_region(context)
                if operators.modal_state.creation.temp_point:
                    return self.start_region(context)
                return {'RUNNING_MODAL'}

            # If not creating, we are trying to start a new region.
            if operators.modal_state.creation.temp_point:
                # A temp_point exists. It's either on the surface, or snapped to an existing point.
                if sid_p <= 1:  # Clicked on surface (-1) or a root point (1)
                    return self.start_region(context)

            return {'RUNNING_MODAL'}
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            if operators.modal_state.extrude_mode:
                return self.end_extrusion(context)
            elif operators.modal_state.scaling_mode:
                return self.end_scaling(context)
            elif operators.modal_state.rotation_mode:
                return self.end_rotation(context)
            elif operators.modal_state.dragging_point:
                return self.end_drag(context)
            elif operators.modal_state.move_subregion_mode:
                return self.end_move_subregion(context)
            elif operators.modal_state.twist_mode:
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
            return self.abort_action(context)
        elif event.type == 'TAB' and event.value == 'PRESS':
            self.end_modal(context)
            return {'FINISHED'}
        
        return {'RUNNING_MODAL'}

    #-- Modal Action Handlers --#
    def start_modal(self, context, event):
        if operators.modal_state.edit_mode:
            self.end_modal(context)
            return {'FINISHED'}
        if context.area.type != 'VIEW_3D':
            return {'CANCELLED'}
        selected = context.active_object
        if not utils.tubegroom_object(selected):
            self.report({'ERROR'}, "Select TubeGroom object")
            return {'CANCELLED'}
        geom.ensure_attrs(selected.data)
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        from . import interpolation
        interpolation.rebuild_regions(selected)
        surf_obj = utils.get_surface_obj(context)
        if not surf_obj or surf_obj.type != 'MESH':
            self.report({'ERROR'}, "Select target object")
            return {'CANCELLED'}

        areas = context.window.screen.areas if context.window else bpy.context.screen.areas
        for area in areas:
            if area.type != 'VIEW_3D':
                continue
            for space in area.spaces:
                if space.type != 'VIEW_3D':
                    continue
                shade = space.shading
                color_type_prop = shade.bl_rna.properties.get('color_type')
                allowed_types = {item.identifier for item in color_type_prop.enum_items} if color_type_prop else set()
                if 'VERTEX' in allowed_types:
                    shade.color_type = 'VERTEX'
        ca = selected.data.color_attributes
        if 'region_color' in ca:
            ca.active_color = ca['region_color']
        if context.mode == 'OBJECT':
            bpy.ops.object.select_all(action='DESELECT')
        operators.clear_modal_state()
        operators.modal_state.edit_mode = True
        operators.modal_state.creation.current_region_color_index = geom.TubeGroom.next_region_id
        operators.save_state()
      
        from . import drawing as draw_regions
        from .drawing import draw_handlers

        draw_handlers.regions_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_regions,
            (context,),
            'WINDOW',
            'POST_VIEW'
        )
        draw_handlers.edges_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_edges,
            (context,),
            'WINDOW',
            'POST_VIEW'
        )
        draw_handlers.points_handler = bpy.types.SpaceView3D.draw_handler_add(
            draw_regions.draw_points,
            (context,),
            'WINDOW',
            'POST_VIEW'
        )

        context.window_manager.modal_handler_add(self)
        self.report({'INFO'}, "Edit mode")
        return {'RUNNING_MODAL'}

    def start_region(self, context):
        point_count = operators.add_point(context)
        self.report({'INFO'}, f"Point {point_count} added")
        return {'RUNNING_MODAL'}
    def end_region(self, context):
        region_id = operators.end_region(context)
        if region_id is None:
            self.report({'WARNING'}, 'Need 3+ points')
            return {'RUNNING_MODAL'}
        self.report({'INFO'}, f"Region {region_id} closed")
        return {'RUNNING_MODAL'}

    def start_drag_point(self, context, event, region_id, subregion_id, point_index):
        if operators.start_drag(context, event, region_id, subregion_id, point_index):
            self.report({'INFO'}, "Dragging point")
        return {'RUNNING_MODAL'}
    def end_drag(self, context):
        operators.end_drag(context)
        self.report({'INFO'}, "Moved")
        return {'RUNNING_MODAL'}

    def start_extrusion(self, context, event):
        ok, region_id, subregion_id, start_mouse_3d = operators.start_extrusion(context, event)
        
        if ok:
            operators.modal_state.drag.start_mouse_3d = start_mouse_3d
            operators.modal_state.selection.region_id = region_id
            
            self.report({'INFO'}, "Extrusion started")
        else:
            self.report({'WARNING'}, "No subregion to extrude")
        
        return {'RUNNING_MODAL'}

    def end_extrusion(self, context):
        region_id = operators.modal_state.selection.region_id
        operators.end_extrusion(context, region_id)
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
        region = geom.TubeGroom.regions.get(region_id)
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
    
    def start_twist(self, context, event, region_id, subregion_id, hierarchical=False):
        operators.start_twist(context, region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Twisting hierarchy' if hierarchical else 'Twisting'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_twist(self, context):
        operators.end_twist(context)
        self.report({'INFO'}, 'Twist done')
        return {'RUNNING_MODAL'}

    def start_scaling(self, context, event, region_id, subregion_id, hierarchical=False):
        operators.start_scaling(region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Scaling hierarchy' if hierarchical else 'Scaling'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_scaling(self, context):
        region_id = operators.modal_state.selection.region_id
        operators.end_scaling(context, region_id)
        self.report({'INFO'}, "Scaling done")
        context.area.tag_redraw()
        return {'RUNNING_MODAL'}
    
    def start_rotation(self, context, event, region_id, subregion_id, hierarchical=False):
        operators.start_rotation(region_id, subregion_id, hierarchical=hierarchical)
        msg = 'Rotating hierarchy' if hierarchical else 'Rotating'
        self.report({'INFO'}, msg)
        return {'RUNNING_MODAL'}
    def end_rotation(self, context):
        region_id = operators.modal_state.selection.region_id
        subregion_id = operators.modal_state.selection.subregion_id
        operators.end_rotation(context, region_id, subregion_id)
        self.report({'INFO'}, "Rotation done")
        context.area.tag_redraw()
        return {'RUNNING_MODAL'}

    def do_undo(self, context):
        return self.history_action(context, operators.undo_state, "Undone", "Nothing to undo")
    def do_redo(self, context):
        return self.history_action(context, operators.redo_state, "Redone", "Nothing to redo")
    def history_action(self, context, action_func, success_msg, fail_msg):
        if action_func():
            geom.update_mesh_date(context)
            from . import interpolation
            interpolation.update_interpolation(context, None, update_topology=True)
            from . import drawing
            drawing.clear_cache()
            context.area.tag_redraw()
            self.report({'INFO'}, success_msg)
        else:
            self.report({'INFO'}, fail_msg)
        return {'RUNNING_MODAL'}
    
    def abort_action(self, context):
        result = operators.cancel_active_operation(context)
        if result:
            self.report({'INFO'}, result)
            context.area.tag_redraw()
        return {'RUNNING_MODAL'}

    def end_modal(self, context):
        from .drawing import remove_draw_handlers
        remove_draw_handlers()
        operators.clear_modal_state()
        from . import drawing as draw_regions
        draw_regions.clear_cache()
        operators.modal_state.edit_mode = False
        if context.area:
            context.area.tag_redraw()
        total_points = len(operators.modal_state.creation.current_region_points)
        for region in geom.TubeGroom.regions.values():
            for sub in region.subregions.values():
                total_points += len(sub.points)
        self.report({'INFO'}, f"{len(geom.TubeGroom.regions)} regions, {total_points} points")
        return {'FINISHED'}
class TUBEGROOM_OT_clear_points(bpy.types.Operator):
    bl_idname = "strand.clear_points"
    bl_label = "Clear All Regions"
    bl_description = "Delete all TubeGroom regions and clear all associated data"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        if not operators.modal_state.edit_mode:
            from . import drawing
            drawing.remove_draw_handlers()
        if context.mode == 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')
        operators.save_state()
        operators.reset_all_data(keep_edit_mode=operators.modal_state.edit_mode)

        for obj in bpy.data.objects:
            if not utils.tubegroom_object(obj):
                continue

            if obj.type == 'MESH' and obj.data:
                mesh = obj.data
                materials = list(mesh.materials)
                mesh.clear_geometry()
                for mat in materials:
                    mesh.materials.append(mat)
                geom.ensure_attrs(mesh)

            base_name = utils.get_base_name(obj)
            if base_name:
                cobj = bpy.data.objects.get(f"CRV_{base_name}")
                if cobj:
                    data_ref = cobj.data
                    bpy.data.objects.remove(cobj, do_unlink=True)
                    bpy.data.hair_curves.remove(data_ref)
        if context.area:
            context.area.tag_redraw()
        self.report({'INFO'}, "Cleared")
        return {'FINISHED'}

# Interpolation Panel
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
        box.prop(scene, 'strand_curves_enabled', text='Enable Curves Display')
        # Sub-layout enabled when the main toggle is on
        sub_box = box.box()
        sub_box.active = scene.strand_curves_enabled
        sub_box.prop(scene, 'strand_interpolation_enabled', text='Live Updates')
        sub_box.prop(scene, 'strand_density_factor', text='Density')
        # Shape settings collapsible section
        row = sub_box.row()
        row.prop(scene, "strand_show_shape_settings", text=" Shape",
                 icon='TRIA_RIGHT' if not scene.strand_show_shape_settings else 'TRIA_DOWN',
                 emboss=False)
        if scene.strand_show_shape_settings:
            shape_box = sub_box.box()
            # These controls mimic the old hair display settings for the curves object
            shape_box.label(text="Viewport Display")
            row = shape_box.row()
            row.prop(scene.render, "hair_type", expand=True)
            shape_box.prop(scene.render, "hair_subdiv", text="Additional Subdivisions")

# Controls Panel
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
            "Shift+Click: Move (edge)(point) / hierarchy (body)",
            "Shift+Alt+Click: Move column (point)",
            "Alt+Click: Scale (edge) / hierarchy (body)",
            "Shift+Alt+Click: Rotate (edge) / hierarchy (body)",
            "Ctrl+Alt+Click: Twist (edge) / hierarchy (body)",
            "Right Click: Delete",
            "Ctrl+Z: Undo",
            "Ctrl+Shift+Z: Redo",
            "Tab: Exit",
        ]
        for text in controls:
            command, _, description = text.partition(":")
            row = col.row(align=True)
            split = row.split(factor=0.45, align=True)
            split.label(text=command)
            split.label(text=description.strip())

# Classes
classes = [
    TUBEGROOM_PT_main_panel,
    TUBEGROOM_PT_mesh_panel,
    TUBEGROOM_PT_interpolation_panel,
    TUBEGROOM_PT_controls_panel,
    TUBEGROOM_OT_create_tubegroom,
    TUBEGROOM_OT_edit_mode,
    TUBEGROOM_OT_mesh_to_tubegroom,
    TUBEGROOM_OT_clear_points,
]
