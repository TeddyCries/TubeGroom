import bpy
import math
from mathutils import Vector
from . import core as shared_data, utils

# Mesh and attribute management
def ensure_region_attributes(mesh):
    attrs = mesh.attributes
    if 'region_id' not in attrs:
        attrs.new(name='region_id', type='INT', domain='POINT')
    if 'subregion_id' not in attrs:
        attrs.new(name='subregion_id', type='INT', domain='POINT')
    if 'region_color' not in mesh.color_attributes:
        mesh.color_attributes.new(name='region_color', type='FLOAT_COLOR', domain='POINT')

def get_or_create_tubegroom_object(context, target_obj=None, allow_create=True):
    active_obj = context.active_object
    if active_obj and utils.is_tubegroom_object(active_obj):
        ensure_region_attributes(active_obj.data)
        shared_data.set_history_owner(active_obj, context.scene)
        return active_obj
    
    if not target_obj:
        target_obj = context.scene.strand_raycast_target
    if not target_obj or target_obj.type != 'MESH' or utils.is_tubegroom_object(target_obj):
        return None
    
    geo_name = f"GEO_TubeGroom_{target_obj.name}"
    obj = bpy.data.objects.get(geo_name)
    if obj:
        ensure_region_attributes(obj.data)
        shared_data.set_history_owner(obj, context.scene)
        return obj
    
    if not allow_create:
        return None
    
    mesh = bpy.data.meshes.new(geo_name)
    obj = bpy.data.objects.new(geo_name, mesh)
    context.scene.collection.objects.link(obj)
    ensure_region_attributes(mesh)
    shared_data.set_history_owner(obj, context.scene)
    
    prev_active = context.view_layer.objects.active
    prev_selected = list(context.selected_objects)
    for o in prev_selected:
        o.select_set(False)
    obj.select_set(True)
    context.view_layer.objects.active = obj
    bpy.ops.object.shade_auto_smooth(use_auto_smooth=True, angle=math.radians(30))
    for o in prev_selected:
        o.select_set(True)
    context.view_layer.objects.active = prev_active
    
    return obj

def update_region_mesh(context, _region_id):
    return create_or_update_merged_mesh(context)

def build_region_mesh_data(context, region_id):
    bend, segs = utils.get_view_settings(context)
    region = shared_data.regions.get(region_id)
    if not region or not region.subregions:
        return [], [], [], [], set()
    
    color = utils.get_unique_color(region.color_index)
    color_tuple = (color[0], color[1], color[2], 1.0)
    rid = region.region_id
    
    real_rings_pos = [([p.position for p in sub.points], (rid, sid)) for sid, sub in sorted(region.subregions.items())]
    rings = utils.calculate_interpolated_rings(real_rings_pos, bend, segs)
    
    vertices, faces, vert_meta, vert_colors = [], [], [], []
    vertical_edges = set()
    vertex_offset = 0
    prev_n, prev_start = None, None
    last_idx = len(rings) - 1
    
    for idx, (ring_pos, meta) in enumerate(rings):
        n = len(ring_pos)
        if n < 3:
            continue
        
        start = vertex_offset
        vertices.extend((p.x, p.y, p.z) for p in ring_pos)
        vert_meta.extend([meta] * n)
        vert_colors.extend([color_tuple] * n)
        
        if prev_start is not None and prev_n == n:
            faces.extend([prev_start + j, prev_start + (j + 1) % n, start + (j + 1) % n, start + j] for j in range(n))
            vertical_edges.update((prev_start + j, start + j) for j in range(n))
        
        if idx == last_idx and meta[1] > 0:
            faces.append(list(range(start, start + n)))
        
        vertex_offset += n
        prev_start, prev_n = start, n
    
    return vertices, faces, vert_meta, vert_colors, vertical_edges

def build_merged_mesh_data(context):
    all_verts, all_faces, all_vert_meta, all_vert_colors = [], [], [], []
    all_vertical_edges = set()
    vertex_offset = 0
    
    for region_id in sorted(shared_data.regions.keys()):
        verts, faces, meta, colors, v_edges = build_region_mesh_data(context, region_id)
        if not verts:
            continue
        
        all_verts.extend(verts)
        all_vert_meta.extend(meta)
        all_vert_colors.extend(colors)
        
        if vertex_offset:
            all_faces.extend([v_idx + vertex_offset for v_idx in face] for face in faces)
            all_vertical_edges.update((a + vertex_offset, b + vertex_offset) for a, b in v_edges)
        else:
            all_faces.extend(faces)
            all_vertical_edges.update(v_edges)
        
        vertex_offset += len(verts)
    
    return all_verts, all_faces, all_vert_meta, all_vert_colors, all_vertical_edges

def create_or_update_merged_mesh(context, allow_rebuild_from_mesh=False, allow_create=True):
    obj = get_or_create_tubegroom_object(context, allow_create=allow_create)
    if not obj:
        if not allow_create and shared_data.regions and not any(utils.is_tubegroom_object(o) for o in bpy.data.objects):
            shared_data.reset_all_data(clear_history=True)
            from . import drawing
            drawing.clear_cache()
        return False
    
    if allow_rebuild_from_mesh and not shared_data.regions:
        from . import interpolation
        interpolation.rebuild_regions_from_merged_mesh(obj)
    
    mesh = obj.data
    verts, faces, vert_meta, vert_colors, vertical_edges = build_merged_mesh_data(context)
    matrix_inv = obj.matrix_world.inverted()
    local_verts = [matrix_inv @ Vector(v) for v in verts] if verts else []
    topology_changed = len(mesh.vertices) != len(local_verts) or len(mesh.polygons) != len(faces)
    
    if topology_changed:
        materials = list(mesh.materials)
        mesh.clear_geometry()
        if local_verts:
            mesh.from_pydata(local_verts, [], faces)
            for mat in materials:
                mesh.materials.append(mat)
            mesh.update()
            ensure_region_attributes(mesh)
            mesh.color_attributes.active_color = mesh.color_attributes['region_color']
            utils.set_mesh_attribute(mesh, 'region_id', 'INT', 'POINT', [r for r, _ in vert_meta])
            utils.set_mesh_attribute(mesh, 'subregion_id', 'INT', 'POINT', [s for _, s in vert_meta])
            utils.set_mesh_attribute(mesh, 'region_color', 'FLOAT_COLOR', 'POINT', vert_colors)
            if vertical_edges:
                mask = [1 if (e.vertices[0], e.vertices[1]) in vertical_edges or (e.vertices[1], e.vertices[0]) in vertical_edges else 0 for e in mesh.edges]
                utils.set_mesh_attribute(mesh, 'sharp_edge', 'BOOLEAN', 'EDGE', mask)
            if 'sharp_face' in mesh.attributes:
                mesh.attributes.remove(mesh.attributes['sharp_face'])
    elif local_verts:
        mesh.vertices.foreach_set('co', [coord for v in local_verts for coord in v])
        mesh.update()
    
    return True


