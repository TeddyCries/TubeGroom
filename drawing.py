import bpy
import gpu
from mathutils import Vector
from gpu_extras.batch import batch_for_shader
from . import utils

# Draw Handlers and Cache
class DrawHandlers:
    def __init__(self):
        self.points_handler = None
        self.edges_handler = None
        self.regions_handler = None
        self.cache = {'regions': {}, 'edges': {}, 'points': {}}

draw_handlers = DrawHandlers()

def remove_draw_handlers():
        if draw_handlers.points_handler:
            bpy.types.SpaceView3D.draw_handler_remove(draw_handlers.points_handler, 'WINDOW')
            draw_handlers.points_handler = None
        if draw_handlers.edges_handler:
            bpy.types.SpaceView3D.draw_handler_remove(draw_handlers.edges_handler, 'WINDOW')
            draw_handlers.edges_handler = None
        if draw_handlers.regions_handler:
            bpy.types.SpaceView3D.draw_handler_remove(draw_handlers.regions_handler, 'WINDOW')
            draw_handlers.regions_handler = None

        clear_cache()

# Cache Utilities
def clear_cache():
        for section in draw_handlers.cache.values():
            section.clear()
def cache_handler(cache_section, cache_key, key, fn):
    cache_dict = draw_handlers.cache[cache_section]
    entry = cache_dict.get(cache_key)
    if not entry or entry['key'] != key:
        cache_dict[cache_key] = {'key': key, 'positions': fn()}
    return cache_dict[cache_key]['positions']
    
# Draw batch utility
def draw_batch(shader, mode, positions, color=None, size=None, width=None):
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

# Geometry Utils
def triangulate_polygon(positions):
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
    idx_tris = utils.tri_face(uv)
    verts = [(p.x, p.y, p.z) for p in positions]
    return [verts[i] for tri in idx_tris if len(tri) == 3 for i in tri]

# Drawing Functions
def draw_regions(context):
    from . import geometry, operators
    obj = utils.get_surface_obj(context)
    root_offset = getattr(context.scene, 'strand_root_offset', 0.005)
    all_tris, all_colors = [], []

    for region_id, region in geometry.regions.items():
        if len(region.subregions) != 1 or 1 not in region.subregions:
            continue
        root = region.subregions.get(1)
        if not root or len(root.get_positions()) < 3:
            continue
        extra = ('rootfill', round(float(root_offset), 6), 4, getattr(obj, 'name', None))
        key = (root.last_modified, extra)
        entry = draw_handlers.cache['regions'].get(region_id)

        if not entry or entry['key'] != key:
            shrink_positions = utils.snap_points(obj, root.get_positions(), root_offset, 4)
            tris = triangulate_polygon(shrink_positions)
            draw_handlers.cache['regions'][region_id] = {'key': key, 'tris': tris}
        else:
            tris = entry['tris']

        if tris:
            rgba = (*utils.color_id(region.color_index), 0.5)
            all_tris.extend(tris)
            all_colors.extend([rgba]*len(tris))

    if len(operators.modal_state.current_region_points) >= 3 and operators.modal_state.current_region_color_index is not None:
        shrink_positions = utils.snap_points(obj, operators.modal_state.current_region_points, root_offset, 4)
        tris = triangulate_polygon(shrink_positions)
        if tris:
            rgba = utils.color_id(operators.modal_state.current_region_color_index)
            all_tris.extend(tris)
            all_colors.extend([rgba]*len(tris))

    if all_tris:
        gpu.state.blend_set('ALPHA')
        for col, tris_chunk in utils.regions_color(all_tris, all_colors).items():
            draw_batch(gpu.shader.from_builtin('UNIFORM_COLOR'), 'TRIS', tris_chunk, col)
        gpu.state.blend_set('NONE')
def draw_edges(context):
    from . import geometry, operators
    obj = utils.get_surface_obj(context)
    root_offset = getattr(context.scene, 'strand_root_offset', 0.005)
    all_edges = []

    for region in geometry.regions.values():
        for sub_id, sub in region.subregions.items():
            if len(sub.get_positions()) < 3:
                continue
            if sub_id == 1:
                extra = ('root_edges', round(float(root_offset), 6), 4, getattr(obj, 'name', None))
                key = (sub.last_modified, extra)
                shrink_positions = cache_handler(
                    'edges', (region.region_id, sub_id), key,
                    lambda: utils.snap_points(obj, sub.get_positions(), root_offset, 4)
                )
            else:
                shrink_positions = sub.get_positions()
            for i in range(len(shrink_positions)):
                start = shrink_positions[i]
                end = shrink_positions[(i + 1) % len(shrink_positions)]
                all_edges.extend((start, end))

    if all_edges:
        gpu.state.depth_test_set('LESS_EQUAL')
        draw_batch(gpu.shader.from_builtin('UNIFORM_COLOR'), 'LINES', all_edges, (1,1,1,1), width=2.0)
        gpu.state.depth_test_set('NONE')

    if len(operators.modal_state.current_region_points) >= 2:
        current_edges = [p for i,p in enumerate(operators.modal_state.current_region_points[:-1])
                         for p in (operators.modal_state.current_region_points[i], operators.modal_state.current_region_points[i+1])]
        gpu.state.depth_test_set('ALWAYS')
        draw_batch(gpu.shader.from_builtin('UNIFORM_COLOR'), 'LINES', current_edges, (1.0,1.0,0.2,1.0), width=2.0)
        gpu.state.depth_test_set('NONE')

    if operators.modal_state.current_region_points and operators.modal_state.temp_point and not operators.modal_state.dragging_point:
        start = operators.modal_state.current_region_points[-1]
        end = operators.modal_state.current_region_points[0] if operators.modal_state.snap_to_first and len(operators.modal_state.current_region_points)>=3 else operators.modal_state.temp_point
        draw_batch(gpu.shader.from_builtin('UNIFORM_COLOR'), 'LINES', [start, end], (1.0,1.0,0.0,0.8), width=2.0)
def draw_points(context):
    from . import geometry, operators
    obj = utils.get_surface_obj(context)
    root_offset = getattr(context.scene, 'strand_root_offset', 0.005)
    points_normal, points_dragged, points_creating, points_creating_dragged, points_highlight = [], [], [], [], []

    for region_id, region in geometry.regions.items():
        for sub_id, sub in region.subregions.items():
            positions = sub.get_positions()
            if not positions:
                continue
            if sub_id == 1:
                extra = ('root_points', round(float(root_offset), 6), 0, getattr(obj, 'name', None))
                key = (sub.last_modified, extra)
                display_positions = cache_handler(
                    'points', (region.region_id, sub_id), key,
                    lambda: utils.snap_points(obj, positions, root_offset, 0)
                )
            else:
                display_positions = positions

            for i, pos in enumerate(display_positions):
                if (operators.modal_state.dragging_point and
                    operators.modal_state.selected_region_id == region_id and
                    operators.modal_state.selected_subregion_id == sub_id and
                    operators.modal_state.selected_point_index == i):
                    points_dragged.append(pos)
                else:
                    if hasattr(operators.modal_state, 'highlight_region_id') and operators.modal_state.highlight_region_id == region_id and operators.modal_state.highlight_subregion_id == sub_id:
                        points_highlight.append(pos)
                    else:
                        points_normal.append(pos)

    if operators.modal_state.current_region_points:
        preview = utils.snap_points(obj, operators.modal_state.current_region_points, root_offset, 0)
        for j, pos in enumerate(preview):
            if operators.modal_state.dragging_point and operators.modal_state.selected_region_id == -1 and operators.modal_state.selected_point_index == j:
                points_creating_dragged.append(pos)
            else:
                points_creating.append(pos)

    gpu.state.depth_test_set('LESS_EQUAL')
    draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', points_normal, (0.2,0.8,1.0,1.0), size=10 if operators.modal_state.edit_mode else 8)
    draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', points_highlight, (1,1,1,1), size=10 if operators.modal_state.edit_mode else 8)
    draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', points_creating, (0.2,0.8,1.0,1.0), size=8)
    gpu.state.depth_test_set('ALWAYS')
    draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', points_dragged, (1,1,1,1), size=16 if operators.modal_state.snap_to_nearest else 14)
    gpu.state.depth_test_set('LESS_EQUAL')

    if points_creating_dragged:
        color = (1,1,0,1) if operators.modal_state.snap_to_nearest else (1,0,1,1)
        size = 16 if operators.modal_state.snap_to_nearest else 14
        gpu.state.depth_test_set('ALWAYS')
        draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', points_creating_dragged, color, size=size)
        gpu.state.depth_test_set('LESS_EQUAL')

    if operators.modal_state.temp_point and not operators.modal_state.snap_to_first:
        color = (1,1,0,1) if operators.modal_state.snap_to_nearest else (0.2,0.8,1.0,0.7)
        size = 12 if operators.modal_state.snap_to_nearest else 8
        gpu.state.depth_test_set('ALWAYS')
        draw_batch(gpu.shader.from_builtin('POINT_UNIFORM_COLOR'), 'POINTS', [operators.modal_state.temp_point], color, size=size)
        gpu.state.depth_test_set('LESS_EQUAL')

    gpu.state.depth_test_set('NONE')
