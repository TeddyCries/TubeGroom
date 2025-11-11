import bpy
from mathutils import Vector
import time
from . import utils, geometry

# Data classes
class Point:
    """Represent a single vertex in a subregion."""
    def __init__(self, position, region_id, subregion_id):
        self.position = position
        self.region_id = region_id
        self.subregion_id = subregion_id

class Subregion:
    """Represent a cross-section of a guide tube."""
    def __init__(self, region_id, subregion_id):
        self.region_id = region_id
        self.subregion_id = subregion_id
        self.points = []
        self.extrude_height = 0.0
        self.last_modified = time.time()
    def add_point(self, position):
        point = Point(position, self.region_id, self.subregion_id)
        self.points.append(point)
        self.last_modified = time.time()
        return point
    def get_positions(self):
        return [p.position for p in self.points]
    def touch(self):
        self.last_modified = time.time()

class Region:
    """Represent a full guide tube with subregions."""
    def __init__(self, region_id, color_index=0):
        self.region_id = region_id
        self.subregions = {}
        self.next_subregion_id = 1
        self.color_index = color_index
        self.last_modified = time.time()
    def create_subregion(self):
        subregion = Subregion(self.region_id, self.next_subregion_id)
        self.subregions[self.next_subregion_id] = subregion
        self.next_subregion_id += 1
        return subregion
    def renumber_subregions(self):
        if not self.subregions:
            self.next_subregion_id = 1
            return
        new_subregions = {}
        new_sid = 1
        for old_sid in sorted(self.subregions.keys()):
            sub = self.subregions[old_sid]
            sub.subregion_id = new_sid
            for p in sub.points:
                p.subregion_id = new_sid
            new_subregions[new_sid] = sub
            new_sid += 1
        self.subregions = new_subregions
        self.next_subregion_id = new_sid
    def touch(self):
        self.last_modified = time.time()

class Face:
    """Represent a face with cached orientations and positions."""
    def __init__(self, face_id, center, normal, tangent, bitangent, positions):
        self.face_id = face_id
        self.subregion_id = face_id
        self.center = Vector(center)
        n = Vector(normal)
        self.normal = n.normalized() if n.length else Vector((0, 0, 1))
        t = Vector(tangent)
        self.tangent = t.normalized() if t.length else Vector((1, 0, 0))
        b = Vector(bitangent)
        self.bitangent = b.normalized() if b.length else self.normal.cross(self.tangent)
        self.region_positions = [Vector(p) for p in positions]
        self.reference_vertices2d = []
        self.tri_indices = []
    def uv_to_world_abs(self, u_abs, v_abs):
        if not self.reference_vertices2d or not self.tri_indices:
            return self.center
        p2d = (u_abs, v_abs)
        for (i, j, k) in self.tri_indices:
            a = self.reference_vertices2d[i]
            b = self.reference_vertices2d[j]
            c = self.reference_vertices2d[k]
            w1, w2, w3, inside = utils.barycentric(p2d, a, b, c)
            if inside:
                a_val = self.region_positions[i]
                b_val = self.region_positions[j]
                c_val = self.region_positions[k]
                return a_val * w1 + b_val * w2 + c_val * w3
        return self.center

# Global state
regions = {}
next_region_id = 1
current_region_points = []
current_region_color_index = 0

temp_point = None
snap_to_first = False
snap_to_existing = None

edit_mode = False
exit_request = False
dragging_point = False
selected_region_id = -1
selected_subregion_id = -1
selected_point_index = -1
drag_start_mouse_3d = None
drag_original_position = None

extrude_mode = False
extrude_preview_subregion_id = None
extrude_normal = None
scaling_mode = False
scaling_pivot_centers = None
scaling_original_positions = None
rotation_mode = False
rotation_pivot_center = None
rotation_original_positions = None
move_subregion_mode = False
move_subregion_region_id = -1
move_subregion_original_positions = None
column_drag_mode = False
column_drag_original_positions = None
column_drag_active_subregion = -1
twist_mode = False
twist_original_positions = None
twist_anchor = None
twist_axis = None

history_stack = []
redo_stack = []
max_history_size = 50
DEFAULT_HISTORY_KEY = "__global__"
_history_registry = {}
_redo_registry = {}
_active_history_key = DEFAULT_HISTORY_KEY

def _clear_history_registry():
    global history_stack, redo_stack, _history_registry, _redo_registry, _active_history_key
    history_stack = []
    redo_stack = []
    _history_registry = {DEFAULT_HISTORY_KEY: history_stack}
    _redo_registry = {DEFAULT_HISTORY_KEY: redo_stack}
    _active_history_key = DEFAULT_HISTORY_KEY

_clear_history_registry()

points_handler = None            # Viewport draw handler for points.
edges_handler = None            # Viewport draw handler for edges.
regions_handler = None           # Viewport draw handler for filled regions.
custom_icons = None            # Preview collection for custom icons.
tubegroom_data = None           # Cached interpolation system data.

def _normalize_history_key(owner, scene=None):
    scene_key = None
    if scene is not None:
        scene_key = getattr(scene, "name_full", None) or getattr(scene, "name", None)
    if scene_key is None and hasattr(owner, "users_scene"):
        scenes = owner.users_scene
        if scenes:
            scene_candidate = scenes[0]
            scene_key = getattr(scene_candidate, "name_full", None) or getattr(scene_candidate, "name", None)
    if scene_key is None:
        import bpy
        ctx_scene = getattr(bpy.context, "scene", None)
        if ctx_scene:
            scene_key = getattr(ctx_scene, "name_full", None) or getattr(ctx_scene, "name", None)
    object_key = None
    if hasattr(owner, "name"):
        base = utils.get_base_name_from_obj(owner)
        name_part = base or getattr(owner, "name", None)
        pointer_fn = getattr(owner, "as_pointer", None)
        pointer_suffix = ""
        if callable(pointer_fn):
            pointer_suffix = f":{pointer_fn():x}"
        if name_part:
            object_key = f"{name_part}{pointer_suffix}"
    elif isinstance(owner, str) and owner:
        object_key = owner
    if not object_key:
        object_key = DEFAULT_HISTORY_KEY
    if scene_key:
        return f"{scene_key}::{object_key}"
    return object_key

def set_history_owner(owner, scene=None):
    global history_stack, redo_stack, _active_history_key
    key = _normalize_history_key(owner, scene)
    _active_history_key = key
    stack = _history_registry.get(key)
    if stack is None:
        stack = []
        _history_registry[key] = stack
    redo = _redo_registry.get(key)
    if redo is None:
        redo = []
        _redo_registry[key] = redo
    history_stack = stack
    redo_stack = redo

#-- State Management (Undo/Redo) --#
def snapshot_state():
    """Captures the current state of all regions and points for the undo system."""
    return {
        'regions': {region_id: {
            'region_id': region.region_id,
            'color_index': region.color_index,
            'next_subregion_id': region.next_subregion_id,
            'subregions': {sub_id: {
                'region_id': sub.region_id,
                'subregion_id': sub.subregion_id,
                'last_modified': sub.last_modified,
                'extrude_height': sub.extrude_height,
                'points': [{'position': p.position.copy()} for p in sub.points]
            } for sub_id, sub in region.subregions.items()}
        } for region_id, region in regions.items()},
        'current_region_points': current_region_points.copy(),
        'next_region_id': next_region_id,
        'current_region_color_index': current_region_color_index
    }

def save_state():
    """Pushes the current state onto the undo history stack and clears the redo stack."""
    state = snapshot_state()
    history_stack.append(state)
    # Any new action invalidates redo history
    redo_stack.clear()
    if len(history_stack) > max_history_size:
        history_stack.pop(0)

def apply_state(state):
    """Restores the application to a previously saved state."""
    global regions, current_region_points, next_region_id, current_region_color_index, tubegroom_data
    regions.clear()
    for region_id, region_data in state['regions'].items():
        region = Region(region_data['region_id'], region_data['color_index'])
        region.next_subregion_id = region_data['next_subregion_id']
        for sub_id, sub_data in region_data['subregions'].items():
            subregion = Subregion(sub_data['region_id'], sub_data['subregion_id'])
            subregion.extrude_height = sub_data['extrude_height']
            subregion.last_modified = sub_data.get('last_modified', time.time())
            for p_data in sub_data['points']:
                subregion.add_point(p_data['position'])
            region.subregions[sub_id] = subregion
        regions[region_id] = region
        region.touch()
    current_region_points = state['current_region_points']
    next_region_id = state['next_region_id']
    current_region_color_index = state['current_region_color_index']
    clear_modal_state(keep_edit_mode=True)
    tubegroom_data = None
    geometry.create_or_update_merged_mesh(bpy.context)

def restore_state():
    """Undo: move current state to redo, pop last undo state and apply it."""
    if not history_stack or not edit_mode:
        return False
    # Save current for redo
    current = snapshot_state()
    redo_stack.append(current)
    # Apply previous from undo stack
    state = history_stack.pop()
    apply_state(state)
    return True

def redo_state():
    """Redo: move current state to undo, pop last redo state and apply it."""
    if not redo_stack or not edit_mode:
        return False
    # Save current for undo
    current = snapshot_state()
    history_stack.append(current)
    if len(history_stack) > max_history_size:
        history_stack.pop(0)
    # Apply redo state
    state = redo_stack.pop()
    apply_state(state)
    return True

def clear_modal_state(keep_edit_mode=None):
    global dragging_point, move_subregion_mode, extrude_mode, scaling_mode, rotation_mode, twist_mode
    global selected_region_id, selected_subregion_id, selected_point_index
    global drag_start_mouse_3d, drag_original_position
    global extrude_preview_subregion_id, extrude_normal
    global scaling_pivot_centers, scaling_original_positions
    global rotation_pivot_center, rotation_original_positions
    global twist_original_positions, twist_anchor, twist_axis
    global column_drag_mode, column_drag_original_positions, column_drag_active_subregion
    global move_subregion_region_id, move_subregion_original_positions
    global edit_mode, exit_request
    dragging_point = False
    move_subregion_mode = False
    extrude_mode = False
    scaling_mode = False
    rotation_mode = False
    twist_mode = False
    edit_mode = edit_mode if keep_edit_mode is None else keep_edit_mode
    exit_request = False
    selected_region_id = -1
    selected_subregion_id = -1
    selected_point_index = -1
    drag_start_mouse_3d = None
    drag_original_position = None
    move_subregion_region_id = -1
    move_subregion_original_positions = None
    extrude_preview_subregion_id = None
    extrude_normal = None
    scaling_pivot_centers = None
    scaling_original_positions = None
    rotation_pivot_center = None
    rotation_original_positions = None
    twist_original_positions = None
    twist_anchor = None
    twist_axis = None
    column_drag_mode = False
    column_drag_original_positions = None
    column_drag_active_subregion = -1


def reset_all_data(keep_edit_mode=False, clear_history=False):
    global regions, current_region_points, temp_point, snap_to_first, snap_to_existing
    global edit_mode, exit_request, dragging_point, selected_region_id, selected_subregion_id, selected_point_index
    global extrude_mode, scaling_mode, rotation_mode, next_region_id, current_region_color_index
    global drag_start_mouse_3d, drag_original_position, points_handler, edges_handler, regions_handler, tubegroom_data
    global twist_mode, twist_original_positions, twist_axis, twist_anchor
    global column_drag_mode, column_drag_original_positions, column_drag_active_subregion
    global move_subregion_mode, move_subregion_region_id, move_subregion_original_positions
    global scaling_pivot_centers, scaling_original_positions, extrude_preview_subregion_id, extrude_normal
    global rotation_pivot_center, rotation_original_positions
    regions.clear()
    current_region_points.clear()
    temp_point = None
    snap_to_first = False
    snap_to_existing = None
    edit_mode = keep_edit_mode
    exit_request = False
    dragging_point = False
    selected_region_id = -1
    selected_subregion_id = -1
    selected_point_index = -1
    extrude_mode = False
    scaling_mode = False
    rotation_mode = False
    move_subregion_mode = False
    move_subregion_region_id = -1
    move_subregion_original_positions = None
    column_drag_mode = False
    column_drag_original_positions = None
    column_drag_active_subregion = -1
    twist_mode = False
    twist_original_positions = None
    twist_axis = None
    twist_anchor = None
    extrude_preview_subregion_id = None
    extrude_normal = None
    scaling_pivot_centers = None
    scaling_original_positions = None
    rotation_pivot_center = None
    rotation_original_positions = None
    if clear_history:
        _clear_history_registry()
    next_region_id = 1
    current_region_color_index = 0
    drag_start_mouse_3d = None
    drag_original_position = None
    if not keep_edit_mode:
        points_handler = None
        edges_handler = None
        regions_handler = None
    tubegroom_data = None
