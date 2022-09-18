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
from compas_tno.algorithms import q_from_qid
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems.problems import plot_svds
import numpy as np

def update_geometry(form, M):

    M.q = q_from_qid(M.q, M.ind, M.Edinv, M.Ei, M.ph)
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    i = 0
    for edge in form.edges_where({'_is_edge': True}):
        form.edge_attribute(edge, 'q', float(M.q[i]))
        i += 1

    i = 0
    for key in form.vertices():
        form.vertex_attribute(key, 'x', M.X[i, 0])
        form.vertex_attribute(key, 'y', M.X[i, 1])
        form.vertex_attribute(key, 'z', M.X[i, 2])
        i += 1

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

    plot_svds(M)

# # # ---------- EXAMPLE 1

# discretisation = [6, 6]
# form = FormDiagram.create_ortho_form(discretisation=discretisation, fix='all')
# find_show_inds(form)

# # # # ---------- EXAMPLE 2

# discretisation = [4, 12]
# form = FormDiagram.create_circular_radial_form(discretisation=discretisation)
# find_show_inds(form)

# # # # ---------- EXAMPLE 3

# discretisation = [6, 6]
# form = FormDiagram.create_cross_form(discretisation=discretisation, fix='corners')
# find_show_inds(form)

# ---------- EXAMPLE 4

# path = '/Users/mricardo/compas_dev/me/inds/three_legs.json'
# mesh = Mesh.from_json(path)
# form = FormDiagram.from_mesh(mesh)

# supports = [0, 10, 24, 6, 7, 2, 5, 18, 22]
# for key in supports:
#     form.vertex_attribute(key, 'is_fixed', True)

# for key in form.vertices():
#     form.vertex_attribute(key, 'pz', -1.0)

# edges = [(17, 18), (22, 23), (5, 13), (0, 19)]
# # edges = [(17, 18), (17, 25), (25, 29), (8, 29)]
# form.assign_inds(edges=edges)

# faces =[
#     [0, 19, 15, 30, 14, 20, 6],
#     [5, 13, 28, 27, 16, 12, 24],
#     [2, 26, 8, 29, 25, 17, 18],
# ]

# for face in faces:
#     id = form.add_face(face)
#     form.face_attribute(id, '_is_loaded', False)

# for u, v in form.edges():
#     if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
#         form.edge_attribute((u, v), '_is_edge', False)

# find_show_inds(form)

# from compas_tno.problems import initialize_loadpath
# from compas_tno.viewers import Viewer

# initialize_loadpath(form)
# M = initialise_form(form)

# path_lp = '/Users/mricardo/compas_dev/me/inds/three_legs_lp.json'
# form.to_json(path_lp)

# update_geometry(form, M)

# view = Viewer(form)
# view.draw_thrust()
# view.show()

# plotter = TNOPlotter(form)
# plotter.draw_form()
# plotter.draw_force()
# plotter.show()

path = '/Users/mricardo/compas_dev/me/inds/three_legs_simple_form.json'
form = FormDiagram.from_json(path)

find_show_inds(form)
