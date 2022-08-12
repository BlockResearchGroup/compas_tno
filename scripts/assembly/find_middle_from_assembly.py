# assembly find middle

from compas_assembly.datastructures import Assembly
from compas.datastructures import Mesh
from compas_view2.app import App

from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import Vector
from compas.geometry import Line

from compas.datastructures import mesh_offset

app = App()

address = '/Users/mricardo/compas_dev/me/blocks/Assembly_Interfaces.json'
address = '/Users/mricardo/compas_dev/me/blocks/crossvault/cross-assembly-interfaces.json'
assembly = Assembly.from_json(address)

div = 0.1
S = Scale.from_factors([div]*3)
assembly.transform(S)

thk = 0.2  # value in mm

for block in assembly.blocks():
    app.add(block)
    pass

for interface in assembly.interfaces():
    app.add(interface)

# lines = []
# for u, v in assembly.graph.edges():
#     pt1 = assembly.node_point(u)
#     pt2 = assembly.node_point(v)
#     line = Line(pt1, pt2)
#     print(line)
#     # lines.append(line)
#     app.add(line)

# mesh = Mesh.from_lines(lines)

# mesh.transform(S)
# print('1', mesh.bounding_box())

# mesh_centroid = mesh.centroid()
# vector = - Vector(*mesh_centroid)
# print('Mesh centroid:', mesh_centroid)
# print('Mesh vector:', vector)

# T = Translation.from_vector(vector)

# mesh.transform(T)
# print('2', mesh.bounding_box())

# mesh.transform(S)
# print('3', mesh.bounding_box())


# mesh_centroid = mesh.centroid()

# vector = Vector(0, 0, -mesh_centroid[2])

# T2 = Translation.from_vector(vector)

# mesh.transform(T2)

# app.add(assembly)
# app.view.camera.target = mesh_centroid
app.show()

# for node in assembly.graph.nodes():
#     mesh.add_vertex(node)
#     centroid = assembly.node_point(node)
#     print(centroid)
#     mesh.vertex_attributes(node, 'xyz', centroid)

# for interface in assembly.interfaces():
#     print(interface)


# for block in assembly.blocks():
#     print(block)
#     centroid = block.centroid
#     mesh.add_vertex()
