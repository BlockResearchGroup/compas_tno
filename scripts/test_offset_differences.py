

import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.shapes import view_shapes
from copy import deepcopy


from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MultiMeshViewer
from compas_viewers.objectviewer import ObjectViewer
from compas.datastructures import Mesh

from compas.geometry import cross_vectors

from compas_tno.datastructures import MeshDos
import math

# ----------------------------------------------------------------------
# ---------------------- TEST OFFSET IN A GEOMETRY ---------------------
# ----------------------------------------------------------------------

def invert_vector(vector):
    return [-1*vector[0], -1*vector[1], -1*vector[2]]


# ----------------------- Cross Vault Data ---------------------------

thk = 0.5
span = 10.0
k = 1.0
n = 1
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

shape = Shape.from_library(data_shape)


# # ----------------------- Arch Data ---------------------------

# incl = 1.0
# H = 01.0*incl
# L = 2.0
# thk = 0.50
# discretisation = 20
# b = 0.2  # Out of plane dimension  of arch
# t = 1.0
# type_structure = 'arch'
# type_formdiagram = 'arch'

# data_shape = {
#     'type': type_structure,
#     'H': H,
#     'L': L,
#     'thk': thk,
#     'discretisation': discretisation,
#     'b': b,
#     't': 0.0,
#     'x0': 0.0
# }

# shape = Shape.from_library(data_shape)


# ----------------------- Dome Data ---------------------------

thk = 0.50
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_spaced_fd'
discretisation = [8, 20]
gradients = True
n = 1

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 0.0
}

shape = Shape.from_library(data_shape)

# ------ CODE

swt = shape.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

intrados = shape.intrados
extrados = shape.extrados

vertices, faces = intrados.to_vertices_and_faces()
intra_meshdos = MeshDos.from_vertices_and_faces(vertices, faces)

vertices, faces = extrados.to_vertices_and_faces()
extra_meshdos = MeshDos.from_vertices_and_faces(vertices, faces)


for n in [0.10]:
    intra_offset = intra_meshdos.offset_mesh(n, 'up')
    extra_offset = extra_meshdos.offset_mesh(n, 'down')

# Settings to visualise

viewer = ObjectViewer()

settings_original = {
    'color': '#0000FF',
    'edges.width': 1,
    'opacity': 0.8,
    'vertices.size': 0,
    'vertices.on': False,
    'edges.on': True,
    'faces.on': True,
    }

settings_offset = {
        'color': '#999999',
        'edges.width': 1,
        'opacity': 0.5,
        'vertices.size': 0,
        'vertices.on': False,
        'edges.on': True,
        'faces.on': True,
        }

# viewer.add(intra_meshdos, name="Intra-Original", settings=settings_original)
viewer.add(extra_meshdos, name="Extra-Original", settings=settings_original)
# viewer.add(intra_offset, name="Intra-Offset-n", settings=settings_offset)
viewer.add(extra_offset, name="Extra-Offset-n", settings=settings_offset)

data_shape['thk'] -= 2* n
shape_offset = Shape.from_library(data_shape)
intra_offset_anal = shape_offset.intrados
extra_offset_anal = shape_offset.extrados

# viewer.add(intra_offset_anal, name="Intra-Offset-anal", settings=settings_offset)
viewer.add(extra_offset_anal, name="Extra-Offset-anal", settings=settings_offset)



# Calculating differences and plotting

print('extra')
max_diff_extrados = 0.0
diff_extrados = {}
for key in extra_meshdos.vertices():
    diff_extrados[key] = extra_offset.vertex_attribute(key, 'z') - extra_offset_anal.vertex_attribute(key, 'z')
    if abs(diff_extrados[key]) > max_diff_extrados:
        max_diff_extrados = abs(diff_extrados[key])
        print(extra_meshdos.vertex_attribute(key, 'z'))
        print(extra_offset.vertex_attribute(key, 'z'))
        print(extra_offset_anal.vertex_attribute(key, 'z'))
print(max_diff_extrados)

print('intra')
max_diff_intrados = 0.0
diff_intrados = {}
for key in intra_meshdos.vertices():
    diff_intrados[key] = intra_offset.vertex_attribute(key, 'z') - intra_offset_anal.vertex_attribute(key, 'z')
    if abs(diff_intrados[key]) > max_diff_intrados:
        max_diff_intrados = abs(diff_intrados[key])
        print(intra_meshdos.vertex_attribute(key, 'z'))
        print(intra_offset.vertex_attribute(key, 'z'))
        print(intra_offset_anal.vertex_attribute(key, 'z'))
print(max_diff_intrados)

from compas_plotters import MeshPlotter
from compas.utilities import i_to_red
plotter = MeshPlotter(extra_meshdos, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(text={key: round(diff_extrados[key], 2) for key in extra_meshdos.vertices()}, facecolor={key: i_to_red(abs(diff_extrados[key])/max_diff_extrados) for key in extra_meshdos.vertices()})
plotter.show()

plotter = MeshPlotter(intra_meshdos, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(text={key: round(diff_intrados[key], 2) for key in intra_meshdos.vertices()}, facecolor={key: i_to_red(abs(diff_intrados[key])/max_diff_intrados) for key in intra_meshdos.vertices()})
plotter.show()


viewer.show()
