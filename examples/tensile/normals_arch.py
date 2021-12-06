import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_plotters import MeshPlotter


span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'arch'  # write the type of form diagram you want and is in the file shape
type_structure = 'arch'
thk = 0.1
discretisation_shape = discretisation

H = 2.0
L = 5

save = False
solutions = {}

objective = 'max_section'  # try 'max'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 'tub']  # in the futture add 'tlb' as variables
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'
tubmax = 0.5

# Create form diagram

path = compas_tno.get('')
address = os.path.join(path, 'form.json')
form = FormDiagram.from_json(address)

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'H': H,
    'L': L,
    'x0': 0,
    'discretisation': discretisation_shape,
    'b': 0.3,
    't': 0.0,
}

vault = Shape.from_library(data_shape)
vault.ro = 20.0

# from compas.datastructures import mesh_weld
# vault.middle = mesh_weld(vault.middle)

from compas_plotters import MeshPlotter
lines_xy = []
lines_xz = []
lines_normal = []

i = 0

for mesh in [vault.intrados, vault.extrados, form]:#, vault.middle]:
    for u, v in mesh.edges():
        ucoord = mesh.vertex_coordinates(u)
        vcoord = mesh.vertex_coordinates(v)
        lines_xy.append(
            {
                'start': ucoord[:2],
                'end': vcoord[:2]
            }
        )

        lines_xz.append(
            {
                'start': [ucoord[0], ucoord[2]],
                'end': [vcoord[0], vcoord[2]]
            }
        )

# for mesh in [vault.middle]:
#     for key in mesh.vertices():
#         coord = mesh.vertex_coordinates(key)
#         if abs(coord[1] - 0.0) < 10e-3:
#             print(i)
#             n = mesh.vertex_normal(key)
#             p1 = [coord[0] + n[0], coord[2] + n[2]]
#             p0 = [coord[0], coord[2]]
#             lines_normal.append(
#                 {
#                     'start': p1,
#                     'end': p0
#                 }
#             )
#             i += 1

for mesh in [vault.middle]:
    for key in mesh.vertices():
        coord = mesh.vertex_coordinates(key)
        if abs(coord[1] - 0.0) < 10e-3:
            print(i)
            n = mesh.vertex_normal(key)
            p1 = [coord[0] + n[0], coord[2] + n[2]]
            p0 = [coord[0], coord[2]]
            lines_normal.append(
                {
                    'start': p1,
                    'end': p0
                }
            )
            i += 1

# plotter = MeshPlotter(vault.intrados)
# plotter.draw_edges()
# plotter.draw_vertices(text={key: str(key) for key in vault.middle.vertices()})
# plotter.show()

print(lines_normal)

plotter = MeshPlotter(vault.intrados)
# plotter.draw_lines(lines_xz)
plotter.draw_lines(lines_normal)
plotter.show()

viewer = Viewer()
viewer.app.add(mesh)
viewer.show()
