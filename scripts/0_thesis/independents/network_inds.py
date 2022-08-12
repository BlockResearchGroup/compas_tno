from numpy import asarray
from numpy.linalg import svd
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt
import numpy as np
from compas.colors import Color

from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas.datastructures import Mesh
from compas_tno.problems import initialise_form
from compas_tno.utilities import slide_diagram
from compas_tno.algorithms import check_independents
import numpy as np

def find_show_inds(form):
    M = initialise_form(form, find_inds=True, printout=True)
    k1 = len(M.ind)

    n = form.number_of_vertices()
    m = len(list(form.edges_where({'_is_edge': True})))
    nb = len(form.fixed())
    ni = n - nb
    k = len(list(form.edges_where({'is_ind': True})))

    print('Number of edges:', m)
    print('Number of vertices:', n)
    print('Number of free vertices:', ni)

    n, m = M.E.shape
    print('E Matrix shape:', M.E.shape)
    print('E rank:', matrix_rank(M.E))
    print('Number of Independents:', k)

    plotter = TNOPlotter(form)
    plotter.draw_form_independents()
    plotter.draw_supports(color=Color.red(), size=10)
    plotter.show()

    _, s, _ = svd(asarray(M.E))
    print('max/min singular vectors E', max(s), min(s), len(s))
    # print(s)

# # # ---------- EXAMPLE 1

discretisation = [6, 6]
form = FormDiagram.create_ortho_form(discretisation=discretisation, fix='all')
find_show_inds(form)

# # # ---------- EXAMPLE 2

discretisation = [4, 12]
form = FormDiagram.create_circular_radial_form(discretisation=discretisation)
find_show_inds(form)

# # # ---------- EXAMPLE 3

discretisation = [6, 6]
form = FormDiagram.create_cross_form(discretisation=discretisation, fix='corners')
find_show_inds(form)

# ---------- EXAMPLE 4

path = '/Users/mricardo/compas_dev/me/inds/three_legs.json'
mesh = Mesh.from_json(path)
form = FormDiagram.from_mesh(mesh)

supports = [0, 10, 24, 6, 7, 2, 5, 18, 22]
for key in supports:
    form.vertex_attribute(key, 'is_fixed', True)

for key in form.vertices():
    form.vertex_attribute(key, 'pz', -1.0)

for u, v in form.edges():
    if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
        form.edge_attribute((u, v), '_is_edge', False)

find_show_inds(form)

# from compas_tno.problems import initialize_loadpath
# from compas_tno.viewers import Viewer

# initialize_loadpath(form)

# view = Viewer(form)
# view.draw_thrust()
# view.show()
