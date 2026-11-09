import argparse
import glob
import math
import os
import re
import struct
import sys

import bpy
from mathutils import Vector


def parse_arguments():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    parser = argparse.ArgumentParser(description="Import a binary learnSPH VTK frame sequence into Blender.")
    parser.add_argument(
        "--pattern",
        default="build-release/res/assignment5_c_obstacle_dam_break_*.vtk",
        help="Glob pattern for the VTK frame sequence.",
    )
    parser.add_argument("--radius", type=float, default=0.006, help="Display radius of one particle.")
    parser.add_argument("--fps", type=int, default=20, help="Playback/render frame rate.")
    parser.add_argument("--speed-max", type=float, default=6.0, help="Speed mapped to the hottest color.")
    parser.add_argument("--start-frame", type=int, default=0, help="First VTK frame index to import.")
    parser.add_argument("--end-frame", type=int, default=-1, help="Last VTK frame index to import; -1 means all.")
    parser.add_argument("--render", action="store_true", help="Render the imported sequence as PNG frames.")
    parser.add_argument(
        "--output",
        default="//sph_animation/",
        help="Output directory for rendered PNG frames.",
    )
    return parser.parse_args(argv)


def read_vtk(path):
    with open(path, "rb") as stream:
        def read_line():
            return stream.readline().decode("ascii").strip()

        assert read_line().startswith("# vtk DataFile")
        read_line()
        assert read_line() == "BINARY"
        assert read_line() == "DATASET UNSTRUCTURED_GRID"

        points_header = read_line().split()
        assert points_header[0] == "POINTS"
        point_count = int(points_header[1])
        positions = struct.unpack(f">{3 * point_count}d", stream.read(24 * point_count))
        stream.readline()

        cells_header = read_line().split()
        assert cells_header[0] == "CELLS"
        cells_size = int(cells_header[2])
        stream.read(4 * cells_size)
        stream.readline()

        cell_types_header = read_line().split()
        assert cell_types_header[0] == "CELL_TYPES"
        stream.read(4 * int(cell_types_header[1]))
        stream.readline()

        point_data_header = read_line().split()
        assert point_data_header[0] == "POINT_DATA"
        field_header = read_line().split()
        assert field_header[0] == "FIELD"

        fields = {}
        for _ in range(int(field_header[2])):
            field = read_line().split()
            name = field[0]
            component_count = int(field[1])
            value_count = int(field[2])
            assert field[3] == "double"
            component_values = struct.unpack(
                f">{component_count * value_count}d",
                stream.read(8 * component_count * value_count),
            )
            stream.readline()
            fields[name] = component_values

    return point_count, positions, fields


def frame_index(path):
    match = re.search(r"_(\d+)\.vtk$", path)
    return int(match.group(1)) if match else -1


def load_sequence(pattern, start_frame, end_frame):
    paths = sorted(glob.glob(pattern), key=frame_index)
    paths = [path for path in paths if "_boundary_" not in os.path.basename(path)]
    if not paths:
        raise RuntimeError(f"No VTK files matched {pattern}")

    if start_frame > 0:
        paths = [path for path in paths if frame_index(path) >= start_frame]
    if end_frame >= 0:
        paths = [path for path in paths if frame_index(path) <= end_frame]
    if not paths:
        raise RuntimeError("The selected frame range is empty.")

    raw_frames = []
    max_points = 0
    bounds_min = Vector((math.inf, math.inf, math.inf))
    bounds_max = Vector((-math.inf, -math.inf, -math.inf))
    for path in paths:
        count, positions, fields = read_vtk(path)
        if "vector" not in fields:
            raise RuntimeError(f"{path} does not contain a velocity vector.")
        velocities = fields["vector"]
        speeds = [
            math.sqrt(
                velocities[3 * i] ** 2
                + velocities[3 * i + 1] ** 2
                + velocities[3 * i + 2] ** 2
            )
            for i in range(count)
        ]
        if count > 0:
            xs = positions[0::3]
            ys = positions[1::3]
            zs = positions[2::3]
            bounds_min.x = min(bounds_min.x, min(xs))
            bounds_min.y = min(bounds_min.y, min(ys))
            bounds_min.z = min(bounds_min.z, min(zs))
            bounds_max.x = max(bounds_max.x, max(xs))
            bounds_max.y = max(bounds_max.y, max(ys))
            bounds_max.z = max(bounds_max.z, max(zs))

        raw_frames.append((count, list(positions), speeds))
        max_points = max(max_points, count)

    if not math.isfinite(bounds_min.x):
        raise RuntimeError("All selected VTK frames contain zero particles.")

    frames = []
    for count, positions, speeds in raw_frames:
        positions.extend([0.0] * (3 * (max_points - count)))
        speeds.extend([0.0] * (max_points - count))
        frames.append((count, positions, speeds))

    return paths, max_points, frames, bounds_min, bounds_max


def look_at(obj, target):
    direction = target - obj.location
    obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()


def make_material(speed_max):
    material = bpy.data.materials.new("SPH velocity")
    material.use_nodes = True
    nodes = material.node_tree.nodes
    links = material.node_tree.links
    nodes.clear()

    output = nodes.new("ShaderNodeOutputMaterial")
    shader = nodes.new("ShaderNodeBsdfPrincipled")
    attribute = nodes.new("ShaderNodeAttribute")
    attribute.attribute_name = "speed"
    attribute.attribute_type = "INSTANCER"
    mapping = nodes.new("ShaderNodeMapRange")
    mapping.inputs["From Min"].default_value = 0.0
    mapping.inputs["From Max"].default_value = speed_max
    mapping.inputs["To Min"].default_value = 0.0
    mapping.inputs["To Max"].default_value = 1.0
    ramp = nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].position = 0.0
    ramp.color_ramp.elements[0].color = (0.02, 0.12, 0.55, 1.0)
    ramp.color_ramp.elements[1].position = 0.45
    ramp.color_ramp.elements[1].color = (0.0, 0.75, 0.9, 1.0)
    third = ramp.color_ramp.elements.new(0.75)
    third.color = (1.0, 0.72, 0.05, 1.0)
    fourth = ramp.color_ramp.elements.new(1.0)
    fourth.color = (0.9, 0.02, 0.02, 1.0)

    links.new(attribute.outputs["Fac"], mapping.inputs["Value"])
    links.new(mapping.outputs["Result"], ramp.inputs["Fac"])
    links.new(ramp.outputs["Color"], shader.inputs["Base Color"])
    links.new(ramp.outputs["Color"], shader.inputs["Emission Color"])
    shader.inputs["Emission Strength"].default_value = 0.35
    shader.inputs["Roughness"].default_value = 0.35
    links.new(shader.outputs["BSDF"], output.inputs["Surface"])
    return material


def make_particle_object(mesh, material):
    obj = bpy.data.objects.new("SPH fluid", mesh)
    bpy.context.collection.objects.link(obj)

    node_group = bpy.data.node_groups.new("SPH particle geometry", "GeometryNodeTree")
    node_group.interface.new_socket(name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    node_group.interface.new_socket(name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")

    nodes = node_group.nodes
    links = node_group.links
    group_input = nodes.new("NodeGroupInput")
    group_output = nodes.new("NodeGroupOutput")
    mesh_to_points = nodes.new("GeometryNodeMeshToPoints")
    radius_attribute = nodes.new("GeometryNodeInputNamedAttribute")
    radius_attribute.data_type = "FLOAT"
    radius_attribute.inputs["Name"].default_value = "radius"
    ico_sphere = nodes.new("GeometryNodeMeshIcoSphere")
    ico_sphere.inputs["Radius"].default_value = 1.0
    ico_sphere.inputs["Subdivisions"].default_value = 1
    instance = nodes.new("GeometryNodeInstanceOnPoints")
    set_material = nodes.new("GeometryNodeSetMaterial")
    set_material.inputs["Material"].default_value = material

    links.new(group_input.outputs["Geometry"], mesh_to_points.inputs["Mesh"])
    links.new(radius_attribute.outputs["Attribute"], mesh_to_points.inputs["Radius"])
    links.new(mesh_to_points.outputs["Points"], instance.inputs["Points"])
    links.new(ico_sphere.outputs["Mesh"], instance.inputs["Instance"])
    links.new(radius_attribute.outputs["Attribute"], instance.inputs["Scale"])
    links.new(instance.outputs["Instances"], set_material.inputs["Geometry"])
    links.new(set_material.outputs["Geometry"], group_output.inputs["Geometry"])

    modifier = obj.modifiers.new("SPH particles", "NODES")
    modifier.node_group = node_group
    return obj


def make_scene(frame_count, fps, bounds_min, bounds_max):
    scene = bpy.context.scene
    scene.frame_start = 0
    scene.frame_end = frame_count - 1
    scene.render.fps = fps
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"

    try:
        scene.render.engine = "BLENDER_EEVEE_NEXT"
    except TypeError:
        scene.render.engine = "BLENDER_EEVEE"

    scene.world.color = (0.008, 0.012, 0.025)
    scene.view_settings.view_transform = "AgX"

    center = (bounds_min + bounds_max) * 0.5
    size = max(bounds_max.x - bounds_min.x, bounds_max.y - bounds_min.y, bounds_max.z - bounds_min.z)

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    bpy.context.collection.objects.link(camera)
    camera.location = center + Vector((1.7 * size, -1.6 * size, 1.05 * size))
    camera_data.lens = 50
    look_at(camera, center)
    scene.camera = camera

    sun_data = bpy.data.lights.new("Sun", "SUN")
    sun_data.energy = 3.0
    sun = bpy.data.objects.new("Sun", sun_data)
    bpy.context.collection.objects.link(sun)
    sun.rotation_euler = (0.55, 0.15, 0.7)

    area_data = bpy.data.lights.new("Fill", "AREA")
    area_data.energy = 450.0
    area_data.shape = "DISK"
    area_data.size = 2.0 * size
    area = bpy.data.objects.new("Fill", area_data)
    bpy.context.collection.objects.link(area)
    area.location = center + Vector((-1.2 * size, -1.0 * size, 1.6 * size))
    look_at(area, center)

    bpy.ops.mesh.primitive_plane_add(size=2.5 * size, location=(center.x, center.y, bounds_min.z))
    ground = bpy.context.object
    ground.name = "Tank floor"
    ground_material = bpy.data.materials.new("Floor")
    ground_material.diffuse_color = (0.025, 0.035, 0.05, 1.0)
    ground.data.materials.append(ground_material)

    return scene


def main():
    args = parse_arguments()
    paths, max_points, frames, bounds_min, bounds_max = load_sequence(
        args.pattern, args.start_frame, args.end_frame
    )

    scene = make_scene(len(frames), args.fps, bounds_min, bounds_max)
    material = make_material(args.speed_max)

    mesh = bpy.data.meshes.new("SPH point cloud")
    mesh.vertices.add(max_points)
    speed_attribute = mesh.attributes.new("speed", "FLOAT", "POINT")
    radius_attribute = mesh.attributes.new("radius", "FLOAT", "POINT")
    obj = make_particle_object(mesh, material)
    obj.data.update()

    def update_frame(current_scene):
        frame = max(0, min(current_scene.frame_current, len(frames) - 1))
        count, positions, speeds = frames[frame]
        radii = [args.radius] * count + [0.0] * (max_points - count)
        mesh.vertices.foreach_set("co", positions)
        speed_attribute.data.foreach_set("value", speeds)
        radius_attribute.data.foreach_set("value", radii)
        mesh.update()

    bpy.app.handlers.frame_change_pre.append(update_frame)
    update_frame(scene)

    if args.render:
        scene.render.filepath = os.path.join(args.output, "sph_")
        scene.render.use_file_extension = True
        bpy.ops.render.render(animation=True)

    print(f"Imported {len(paths)} VTK frames with up to {max_points} particles.")
    print(f"Frame range: {scene.frame_start}..{scene.frame_end} at {args.fps} FPS.")


if __name__ == "__main__":
    main()
