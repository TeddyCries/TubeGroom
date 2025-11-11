import gpu
from mathutils import Vector
from gpu_extras.batch import batch_for_shader
from . import core as shared_data
from . import utils

# --- Cache ---
_CACHE = {
    'root_tris': {},      # region_id -> {'key': key, 'tris': [...]}
    'edges': {},          # (region_id, subregion_id) -> {'key': key, 'positions': [...]}
    'points': {},         # (region_id, subregion_id) -> {'key': key, 'positions': [...]}
}

def clear_cache():
    for section in _CACHE.values():
        section.clear()

def remove_draw_handlers():
    import bpy
    if shared_data.points_handler:
        bpy.types.SpaceView3D.draw_handler_remove(shared_data.points_handler, 'WINDOW')
        shared_data.points_handler = None
    if shared_data.edges_handler:
        bpy.types.SpaceView3D.draw_handler_remove(shared_data.edges_handler, 'WINDOW')
        shared_data.edges_handler = None
    if shared_data.regions_handler:
        bpy.types.SpaceView3D.draw_handler_remove(shared_data.regions_handler, 'WINDOW')
        shared_data.regions_handler = None

# --- Shaders ---
_SHADER_COLOR = gpu.shader.from_builtin('UNIFORM_COLOR')
_SHADER_POINT = gpu.shader.from_builtin('POINT_UNIFORM_COLOR')

def draw_batch(shader, mode, positions, color=None, size=None, width=None):
    """Utility for consistent drawing calls."""
    if not positions:
        return
    batch = batch_for_shader(shader, mode, {"pos": positions})
    if size:
        gpu.state.point_size_set(size)
    if width:
        gpu.state.line_width_set(width)
    shader.bind()
    if color is not None:
        if hasattr(color, '__iter__') and not isinstance(color, str):
            col_seq = tuple(float(c) for c in color)
        else:
            col_seq = (float(color),)
        if len(col_seq) == 3:
            col_seq = (*col_seq, 1.0)
        elif len(col_seq) == 1:
            col_seq = col_seq * 4
        elif len(col_seq) > 4:
            col_seq = col_seq[:4]
        shader.uniform_float("color", col_seq)
    batch.draw(shader)
def cache_get_or_create(cache_dict, cache_key, key, fn):
    """Get from cache or create using fn()."""
    entry = cache_dict.get(cache_key)
    if not entry or entry['key'] != key:
        cache_dict[cache_key] = {'key': key, 'positions': fn()}
    return cache_dict[cache_key]['positions']

def _geometry_key(timestamp, extra=None):
    return (timestamp, extra)

# --- Geometry Utils ---
def triangles_from_positions(positions):
    if not positions or len(positions) < 3:
        return []
    c = sum(positions, Vector()) / len(positions)

    def poly_normal(pts):
        acc = Vector((0.0, 0.0, 0.0))
        m = len(pts)
        for i in range(m):
            a = pts[i] - pts[i-1]
            b = pts[(i+1) % m] - pts[i]
            acc += a.cross(b)
        return acc.normalized() if acc.length > 0 else Vector((0,0,1))

    n = poly_normal(positions)
    e = (positions[1] - positions[0]) if len(positions) > 1 else Vector((1,0,0))
    e = e - n * e.dot(n)
    if e.length <= 1e-9:
        tmp = Vector((1,0,0)) if abs(n.x) < 0.9 else Vector((0,1,0))
        e = tmp - n * tmp.dot(n)
    t = e.normalized() if e.length > 1e-9 else Vector((1,0,0))
    b = n.cross(t).normalized()

    uv = [((p - c).dot(t), (p - c).dot(b)) for p in positions]
    idx_tris = utils.triangulate_face(uv)
    verts = [(p.x, p.y, p.z) for p in positions]
    return [verts[i] for tri in idx_tris if len(tri) == 3 for i in tri]

# --- Draw Handlers ---
def draw_regions(self, context):
    obj = utils.get_surface_object(context)
    all_tris, all_colors = [], []

    for region_id, region in shared_data.regions.items():
        if len(region.subregions) != 1 or 1 not in region.subregions:
            continue
        root = region.subregions.get(1)
        if not root or len(root.get_positions()) < 3:
            continue
        extra = ('rootfill', round(float(shared_data.root_offset), 6), 4, getattr(obj, 'name', None))
        key = _geometry_key(root.last_modified, extra)
        entry = _CACHE['root_tris'].get(region_id)

        if not entry or entry['key'] != key:
            shrink_positions = utils.create_shrinkwrap_points(obj, root.get_positions(), shared_data.root_offset, 4)
            tris = triangles_from_positions(shrink_positions)
            _CACHE['root_tris'][region_id] = {'key': key, 'tris': tris}
        else:
            tris = entry['tris']

        if tris:
            rgba = (*utils.get_unique_color(region.color_index), 0.5)
            all_tris.extend(tris)
            all_colors.extend([rgba]*len(tris))

    if len(shared_data.current_region_points) >= 3:
        shrink_positions = utils.create_shrinkwrap_points(obj, shared_data.current_region_points, shared_data.root_offset, 4)
        tris = triangles_from_positions(shrink_positions)
        if tris:
            rgba = utils.get_unique_color(shared_data.current_region_color_index)
            all_tris.extend(tris)
            all_colors.extend([rgba]*len(tris))

    if all_tris:
        gpu.state.blend_set('ALPHA')
        for col, tris_chunk in utils.group_by_color(all_tris, all_colors).items():
            draw_batch(_SHADER_COLOR, 'TRIS', tris_chunk, col)
        gpu.state.blend_set('NONE')

def draw_edges(self, context):
    obj = utils.get_surface_object(context)
    all_edges = []

    for region in shared_data.regions.values():
        for sub_id, sub in region.subregions.items():
            if len(sub.get_positions()) < 3:
                continue
            if sub_id == 1:
                extra = ('root_edges', round(float(shared_data.root_offset), 6), 4, getattr(obj, 'name', None))
                key = _geometry_key(sub.last_modified, extra)
                shrink_positions = cache_get_or_create(
                    _CACHE['edges'], (region.region_id, sub_id), key,
                    lambda: utils.create_shrinkwrap_points(obj, sub.get_positions(), shared_data.root_offset, 4)
                )
            else:
                shrink_positions = sub.get_positions()
            for i in range(len(shrink_positions)):
                start = shrink_positions[i]
                end = shrink_positions[(i + 1) % len(shrink_positions)]
                all_edges.extend((start, end))

    if all_edges:
        gpu.state.depth_test_set('LESS_EQUAL')
        draw_batch(_SHADER_COLOR, 'LINES', all_edges, (1,1,1,1), width=2.0)
        gpu.state.depth_test_set('NONE')

    if len(shared_data.current_region_points) >= 2:
        current_edges = [p for i,p in enumerate(shared_data.current_region_points[:-1])
                         for p in (shared_data.current_region_points[i], shared_data.current_region_points[i+1])]
        gpu.state.depth_test_set('ALWAYS')
        draw_batch(_SHADER_COLOR, 'LINES', current_edges, (1.0,1.0,0.2,1.0), width=2.0)
        gpu.state.depth_test_set('NONE')

    if shared_data.current_region_points and shared_data.temp_point and not shared_data.dragging_point:
        start = shared_data.current_region_points[-1]
        end = shared_data.current_region_points[0] if shared_data.snap_to_first and len(shared_data.current_region_points)>=3 else shared_data.temp_point
        draw_batch(_SHADER_COLOR, 'LINES', [start, end], (1.0,1.0,0.0,0.8), width=2.0)

def draw_points(self, context):
    obj = utils.get_surface_object(context)
    points_normal, points_dragged, points_creating, points_creating_dragged = [], [], [], []

    for region_id, region in shared_data.regions.items():
        for sub_id, sub in region.subregions.items():
            positions = sub.get_positions()
            if not positions:
                continue
            if sub_id == 1:
                extra = ('root_points', round(float(shared_data.root_offset), 6), 0, getattr(obj, 'name', None))
                key = _geometry_key(sub.last_modified, extra)
                display_positions = cache_get_or_create(
                    _CACHE['points'], (region.region_id, sub_id), key,
                    lambda: utils.create_shrinkwrap_points(obj, positions, shared_data.root_offset, 0)
                )
            else:
                display_positions = positions

            for i, pos in enumerate(display_positions):
                if (shared_data.dragging_point and
                    shared_data.selected_region_id == region_id and
                    shared_data.selected_subregion_id == sub_id and
                    shared_data.selected_point_index == i):
                    points_dragged.append(pos)
                else:
                    points_normal.append(pos)

    if shared_data.current_region_points:
        preview = utils.create_shrinkwrap_points(obj, shared_data.current_region_points, shared_data.root_offset, 0)
        for j, pos in enumerate(preview):
            if shared_data.dragging_point and shared_data.selected_region_id == -1 and shared_data.selected_point_index == j:
                points_creating_dragged.append(pos)
            else:
                points_creating.append(pos)

    gpu.state.depth_test_set('LESS_EQUAL')
    draw_batch(_SHADER_POINT, 'POINTS', points_normal, (0.2,0.8,1.0,1.0), size=10 if shared_data.edit_mode else 8)
    draw_batch(_SHADER_POINT, 'POINTS', points_creating, (0.2,0.8,1.0,1.0), size=8)
    draw_batch(_SHADER_POINT, 'POINTS', points_dragged, (1,1,1,1), size=16 if shared_data.snap_to_existing else 14)

    if points_creating_dragged:
        color = (1,1,0,1) if shared_data.snap_to_existing else (1,0,1,1)
        size = 16 if shared_data.snap_to_existing else 14
        draw_batch(_SHADER_POINT, 'POINTS', points_creating_dragged, color, size=size)

    if shared_data.temp_point and not shared_data.snap_to_first and not shared_data.dragging_point:
        color = (1,1,0,1) if shared_data.snap_to_existing else (0.2,0.8,1.0,0.7)
        size = 12 if shared_data.snap_to_existing else 8
        draw_batch(_SHADER_POINT, 'POINTS', [shared_data.temp_point], color, size=size)

    gpu.state.depth_test_set('NONE')
