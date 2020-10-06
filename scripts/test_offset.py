

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


# ----------------------- Arch Data ---------------------------

incl = 1.0
H = 01.0*incl
L = 2.0
thk = 0.50
discretisation = 20
b = 0.2  # Out of plane dimension  of arch
t = 1.0
type_structure = 'arch'
type_formdiagram = 'arch'

data_shape = {
    'type': type_structure,
    'H': H,
    'L': L,
    'thk': thk,
    'discretisation': discretisation,
    'b': b,
    't': 0.0,
    'x0': 0.0
}

shape = Shape.from_library(data_shape)


# ----------------------- Dome Data ---------------------------

# thk = 0.50
# radius = 5.0
# type_structure = 'dome'
# type_formdiagram = 'radial_spaced_fd'
# discretisation = [50, 20]
# gradients = True
# n = 1

# data_shape = {
#     'type': type_structure,
#     'thk': thk,
#     'discretisation': [discretisation[0]*n, discretisation[1]*n],
#     'center': [5.0, 5.0],
#     'radius': radius,
#     't': 0.0
# }

# shape = Shape.from_library(data_shape)

# ------ CODE

swt = shape.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

intrados = shape.intrados
extrados = shape.extrados

vertices, faces = intrados.to_vertices_and_faces()
intra_offset = MeshDos.from_vertices_and_faces(vertices, faces)
# intra_same_xy = MeshDos.from_vertices_and_faces(vertices, faces)

vertices, faces = extrados.to_vertices_and_faces()
extra_offset = MeshDos.from_vertices_and_faces(vertices, faces)
# extra_same_xy = MeshDos.from_vertices_and_faces(vertices, faces)

meshes = [intrados, extrados]

# viewer = MultiMeshViewer()
# viewer.meshes = meshes
# viewer.show()

for n in [0.05]:
    vertices, faces = intra_offset.to_vertices_and_faces()
    intra_offset = MeshDos.from_vertices_and_faces(vertices, faces)

    vertices, faces = extra_offset.to_vertices_and_faces()
    extra_offset = MeshDos.from_vertices_and_faces(vertices, faces)

    print('INTRADOS')
    new_intrados = []
    for key in intra_offset.vertices():
        normal = (intra_offset.vertex_normal(key))
        x, y, z = intra_offset.vertex_coordinates(key)
        deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
        print(n*1/deviation)
        new_intrados.append(z + n*1/deviation)
        # if not intrados.is_vertex_on_boundary(key):
        #     intra_offset.vertex_attribute(key, 'x', x + normal[0]*n)
        #     intra_offset.vertex_attribute(key, 'y', y + normal[1]*n)
        #     intra_offset.vertex_attribute(key, 'z', z + normal[2]*n)
        # else:
        #     deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
        #     intra_offset.vertex_attribute(key, 'z', z + n*1/deviation)
    i = 0
    for key in intra_offset.vertices():
        intra_offset.vertex_attribute(key, 'z', new_intrados[i])
        i += 1
    meshes.append(intra_offset)

    print('EXTRADOS')
    new_extrados = []
    for key in extra_offset.vertices():
        normal = extra_offset.vertex_normal(key)
        x, y, z = extra_offset.vertex_coordinates(key)
        deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
        print(- n*1/deviation)
        new_extrados.append(z - n*1/deviation)
        # print('x y:', x, y, 'normal:', normal)
        # if not extrados.is_vertex_on_boundary(key):
        #     extra_offset.vertex_attribute(key, 'x', x + normal[0]*n)
        #     extra_offset.vertex_attribute(key, 'y', y + normal[1]*n)
        #     extra_offset.vertex_attribute(key, 'z', z + normal[2]*n)
        # else:  # correction for vertices in the boundary for extrados
        #     deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
        #     extra_offset.vertex_attribute(key, 'z', z - n*1/deviation)
    i = 0
    for key in extra_offset.vertices():
        extra_offset.vertex_attribute(key, 'z', new_extrados[i])
        i += 1
    meshes.append(extra_offset)

viewer = MultiMeshViewer()
viewer.meshes = meshes  # [intra_offset, intrados, extra_offset, extrados]
viewer.show()
