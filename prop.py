import bpy
import math

def get_base_name_safe(context):
    """Helper to safely get base_name from TubeGroom object."""
    from . import geometry, utils
    tg_object = geometry.get_tg_obj(context, allow_create=False)
    return utils.get_base_name(tg_object) if tg_object else None
def update_curves_object_color(curves_obj, live_enabled):
    """Helper to update curves object color and display settings."""
    curves_obj.show_in_front = live_enabled
    curves_obj.color[3] = 0.8 if live_enabled else 1.0
def calculate_log_ranges():
    """Helper to calculate logarithmic ranges for radius conversion."""
    return math.log(0.001), math.log(0.5)

# Root offset property update functions
def get_root_offset(self):
    """Expose root offset in millimeters while storing in meters."""
    return self.strand_root_offset * 1000.0
def set_root_offset(self, value):
    """Store root offset in meters while exposing millimeters in the UI."""
    meters = max(0.0, min(0.1, float(value) / 1000.0 if isinstance(value, (int, float)) else 0.0))
    self.strand_root_offset = meters
def refresh_mesh_system(self, context):
    """Generic update function to redraw the TubeGroom mesh."""
    from . import geometry, operators, utils
    
    if operators.modal_state.edit_mode:
        obj = None
        for o in bpy.data.objects:
            if utils.tubegroom_object(o):
                obj = o
                break
        if not obj:
            return
    else:
        active_obj = context.active_object
        if not active_obj or not utils.tubegroom_object(active_obj) or not active_obj.select_get():
            return
        obj = active_obj
    
    if not geometry.TubeGroom.regions and obj.data.vertices:
        from . import interpolation
        interpolation.rebuild_regions(obj)
    geometry.update_mesh_date(context)
    from . import interpolation
    interpolation.update_interpolation(context, None, update_topology=False)

# Interpolation property update functions
def update_curves_visibility(self, context):
    """Shows/hides and rebuilds the curves object if necessary."""
    from . import interpolation
    enabled = context.scene.strand_curves_enabled
    base_name = get_base_name_safe(context)
    if base_name:
        curve_names = [f"CRV_{base_name}"]
    else:
        curve_names = [name for name in bpy.data.objects.keys() if name.startswith("CRV_TubeGroom_")]
    for name in curve_names:
        obj = bpy.data.objects.get(name)
        if obj:
            obj.hide_set(not enabled)
            obj.hide_viewport = not enabled
    if enabled and base_name:
        system = interpolation.generate_interpolation()
        if system:
            interpolation.build_curves_object(base_name, system)
def update_live_interpolation(self, context):
    """Toggle live interpolation and trigger an immediate rebuild when enabling."""
    if not context.scene.strand_curves_enabled:
        return
    enable = context.scene.strand_interpolation_enabled
    base_name = get_base_name_safe(context)
    if not base_name:
        return
    if enable:
        # Build or update curves immediately when enabling live
        from . import interpolation
        interpolation.update_interpolation(context)
    obj = bpy.data.objects.get(f"CRV_{base_name}")
    if obj:
        update_curves_object_color(obj, enable)
def update_interpolation_params(self, context):
    """Updates the interpolation when the radius or density changes."""
    from . import interpolation, utils

    scene = context.scene

    radius = scene.strand_poisson_radius
    clamped_radius = max(0.001, min(0.5, radius))
    log_min, log_max = calculate_log_ranges()
    log_r = math.log(clamped_radius)
    density = 0.5 if abs(log_max - log_min) < 1e-9 else (log_max - log_r) / (log_max - log_min)
    if abs(scene.strand_density_factor - density) > 1e-6:
        scene.strand_density_factor = density

    if not scene.strand_curves_enabled:
        return

    active_obj = context.active_object
    if not active_obj or not utils.tubegroom_object(active_obj):
        return

    base_name = utils.get_base_name(active_obj)
    if not base_name:
        return

    system = interpolation.generate_interpolation()
    if not system:
        return

    curves_obj = interpolation.build_curves_object(base_name, system)
    if curves_obj:
        live_enabled = scene.strand_interpolation_enabled
        update_curves_object_color(curves_obj, live_enabled)
def update_density(self, context):
    """Updates the radius when the density slider changes."""
    scene = context.scene
    density = max(0.0, min(1.0, scene.strand_density_factor))
    log_min, log_max = calculate_log_ranges()
    log_r = log_max - density * (log_max - log_min)
    new_radius = math.exp(log_r)
    if abs(scene.strand_poisson_radius - new_radius) > 1e-6:
        scene.strand_poisson_radius = max(0.001, min(0.5, new_radius))
def update_hair_settings(self, context):
    """Apply render shape settings to TubeGroom curve objects."""
    for obj in bpy.data.objects:
        if obj.name.startswith("CRV_TubeGroom_") and obj.type == 'CURVES':
            obj.data.display.surface_subdivisions = context.scene.render.hair_subdiv
            obj.data.update_tag()

# Property definitions
SCENE_PROPERTIES = [
    ('strand_root_offset', bpy.props.FloatProperty(name="Root Offset", description="Offset applied to the root subregions", default=0.005, min=0.0, max=0.1, step=0.0001, precision=6, update=refresh_mesh_system)),
    ('strand_root_offset_mm', bpy.props.FloatProperty(
        name="Root Offset (mm)",
        description="Offset applied to the root subregions",
        default=5.0,
        min=0.0,
        max=100.0,
        step=0.1,
        precision=3,
        set=set_root_offset,
        get=get_root_offset
    )),
    ('strand_raycast_target', bpy.props.PointerProperty(name="Target Object", description="Object to use as the surface for placing and projecting guides", type=bpy.types.Object, poll=lambda self, obj: obj.type == 'MESH')),
    ('strand_view_segments', bpy.props.IntProperty(name='Resolution', default=8, min=0, max=8, soft_max=8, description='Insert real loops between subregions (0 disables)', update=refresh_mesh_system, options={'HIDDEN'})),
    ('strand_bend_factor', bpy.props.FloatProperty(name='Smooth', default=1.0, min=0.0, max=1.0, subtype='FACTOR', description='Curve only inserted segments between subregions; original subregions stay fixed (0 disables)', update=refresh_mesh_system, options={'HIDDEN'})),
    ('strand_show_shape_settings', bpy.props.BoolProperty(name="Show Shape Settings", description="Show/hide Cycles hair shape settings", default=True)),
    ('strand_curves_enabled', bpy.props.BoolProperty(name="Enable Curves Display", description="Enable, show, and allow recalculation of the final curves object", default=True, update=update_curves_visibility)),
    ('strand_interpolation_enabled', bpy.props.BoolProperty(name='Enable Live Interpolation', default=False, description='If enabled, interpolation updates while editing', update=update_live_interpolation)),
    ('strand_density_factor', bpy.props.FloatProperty(name="Density", description="Controls the density of interpolated curves. Left is less, right is more.", default=0.5, min=0.0, max=1.0, subtype='FACTOR', update=update_density)),    ('strand_poisson_radius', bpy.props.FloatProperty(name="Radius (UV)", default=0.02, min=1e-5, soft_max=0.5, subtype='DISTANCE', step=1, precision=4, update=update_interpolation_params)),
]
RENDER_PROPERTIES = [
    ('hair_type', bpy.props.EnumProperty(name="Hair Shape Type", items=[('STRAND', "Strand", "Display hair as strands"), ('STRIP', "Strip", "Display hair as ribbons")], default='STRAND', update=update_hair_settings)),
    ('hair_subdiv', bpy.props.IntProperty(name="Additional Subdivisions", description="Viewport subdivisions for the curve display", default=0, min=0, max=3, update=update_hair_settings)),
]
