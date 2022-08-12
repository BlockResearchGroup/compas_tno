from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.shapes import MeshDos
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import extended_dual
from compas_tno.utilities import blocks_from_dual

from compas.geometry import Scale
from compas.geometry import Translation

from compas.datastructures import mesh_dual
from compas.datastructures import Mesh
from compas.geometry import add_vectors, scale_vector

from compas.colors import Color

import compas_tno

from compas_assembly.datastructures import Assembly
from compas_assembly.datastructures import Block

from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Vector


from compas_plotters import Plotter

span = 10.0
delta = 1.0
xspan = yspan = [0.0, span]
xspan_vault = yspan_vault = [- delta, span + delta]
thk = 0.50

# jsonpath = '/Users/mricardo/compas_dev/me/pattern/singular/crossvault/mesh-C2.json'
form = FormDiagram.create_cross_form(discretisation=6)
# form = Mesh.from_json(jsonpath)
# form.parameters = {'type': 'general', 't': 0.0}
# form = FormDiagram.create_fan_form()
vaultfull = Shape.create_crossvault(xy_span=[xspan_vault, yspan_vault])
# vaultfull = Shape.create_pointedcrossvault(xy_span=[xspan_vault, yspan_vault], he=[6.0, 6.0, 6.0, 6.0])

if form.parameters['type'] == 'general':

    bbox = form.bounding_box_xy()
    xmin, ymin = min([m[0] for m in bbox]), min([m[1] for m in bbox])

    if abs(xmin) > 10e-3 or abs(ymin) > 10e-3:
        dx, dy = -xmin, -ymin
        translation = Translation.from_vector([dx, dy, 0.0])
        form.transform(translation)

        bbox = form.bounding_box_xy()
        xmax, ymax = max([m[0] for m in bbox]), max([m[1] for m in bbox])
        print(bbox)
        print(xmax, ymax)

        scale = Scale.from_factors([10.0/xmax, 10.0/ymax, 0.0])
        form.transform(scale)

        bbox = form.bounding_box_xy()
        print(bbox)


apply_envelope_from_shape(form, vaultfull)
apply_selfweight_from_shape(form, vaultfull, normalize=False)

for key in form.vertices():
    target = form.vertex_attribute(key, 'target')
    form.vertex_attribute(key, 'z', target)

# path = compas_tno.get('input-form.json')
# form.to_json(path)

form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM-3.json')

dual = extended_dual(form, Mesh)

# idos, edos = offset_dual(dual, thk)
assembly = blocks_from_dual(dual, thk)

# # remove last block (Not sure why it is duplicated)
# assembly = Assembly()
# lastblock = len(list(assembly_all.blocks())) - 1
# for node in assembly_all.nodes():
#     block = assembly_all.node_block(node)
#     if node < lastblock:
#         assembly.add_block(block)

import compas_tno
file = compas_tno.get('assembly-tno.json')
assembly.to_json(file)
print('Saved Assembly to:', file)

plotter = Plotter()
# artist = plotter.add(dual)

text = {}
for node in assembly.nodes():
    print(node)
    block = assembly.node_block(node)
    centroid = block.centroid()
    x, y, z = centroid
    print(node, centroid)
    text[node] = str(node)
    plotter.axes.text(x, y, str(node))
    # artist = plotter.add(Point(*centroid), text=str(node))
    # artist.draw_labels()

# artist.draw_facelabels()
plotter.show()

view = Viewer(form, show_grid=False)
view.settings['size.edge.max_thickness'] = 5.0
view.draw_thrust(scale_width=False)
# view.draw_shape()
# view.draw_mesh(dual)
# view.draw_mesh(idos)
# view.draw_mesh(edos)

# assembly.to_json()

print(form.parameters)

if form.parameters['type'] == 'cross_fd':
    supports = [91, 101, 111, 81]
    # supports = [31, 37, 43, 25]
elif form.parameters['type'] == 'fan_fd':
    supports = [161, 171, 181, 191]
elif form.parameters['type'] == 'general':
    supports = [267, 279, 255, 243]  # A2
    supports = [190, 203, 216, 177]  # C2
elif form.parameters['type'] is None:
    supports = list(range(81, 97))

print(supports)

# view.draw_assembly(assembly)
for node in assembly.nodes():
    block = assembly.node_block(node)
    if node in supports:
        view.add(block, opacity=0.5, color=Color.red())
    else:
        view.add(block, opacity=0.5)
view.show()
