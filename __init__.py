"""
TubeGroom: tools for creating, editing, and managing tube-based geometry in Blender.
"""

bl_info = {
    "name": "TubeGroom",
    "author": "TeddyCry",
    "version": (1, 0, 1),
    "blender": (4, 5, 0),
    "location": "View3D > Sidebar > TubeGroom",
    "description": "Interactive mesh grooming tools",
    "category": "Mesh",
}

def register():
    import bpy
    import os
    from . import prop, ui

    for name, p in prop.SCENE_PROPERTIES:
        setattr(bpy.types.Scene, name, p)
    for name, p in prop.RENDER_PROPERTIES:
        setattr(bpy.types.RenderSettings, name, p)

    icon_path = os.path.join(os.path.dirname(__file__), "icons", "TG_Main_Icon.png")
    if os.path.exists(icon_path):
        import bpy.utils.previews
        if not hasattr(bpy.utils.previews, 'custom_icons'):
            bpy.utils.previews.custom_icons = bpy.utils.previews.new()
        bpy.utils.previews.custom_icons.load("main_icon", icon_path, 'IMAGE')
    for cls in ui.classes:
        bpy.utils.register_class(cls)

def unregister():
    from . import prop, ui
    try:
        ui.operators.modal_state.edit_mode = False
        from . import drawing
        drawing.remove_draw_handlers()
        ui.operators.modal_state.clear_modal_state()
    except Exception:
        pass

    import bpy.utils.previews
    if hasattr(bpy.utils.previews, 'custom_icons'):
        bpy.utils.previews.remove(bpy.utils.previews.custom_icons)
        delattr(bpy.utils.previews, 'custom_icons')
    for cls in reversed(ui.classes):
        bpy.utils.unregister_class(cls)

    for name, _ in reversed(prop.SCENE_PROPERTIES):
        if hasattr(bpy.types.Scene, name):
            delattr(bpy.types.Scene, name)
    for name, _ in reversed(prop.RENDER_PROPERTIES):
        if hasattr(bpy.types.RenderSettings, name):
            delattr(bpy.types.RenderSettings, name)
