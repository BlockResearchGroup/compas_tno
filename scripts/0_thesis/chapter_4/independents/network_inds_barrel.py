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
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.problems.problems import plot_svds
from compas_tno.problems import initialize_loadpath

from compas.numerical import nullspace as matrix_nullspace
from numpy import array

from compas.colors import Color

from numpy.random import random_sample
from compas_tno.viewers import Viewer
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
    print('Number of Inext.mechanisms:', k-(m-n))

    plotter = TNOPlotter(form)
    plotter.draw_form_independents()
    plotter.draw_supports(color=Color.red(), size=10)
    plotter.show()

    _, s, _ = svd(asarray(M.E))
    print('max/min singular vectors E', max(s), min(s), len(s))
    # print(s)

    plot_svds(M)

def find_nullspaces(form, M):

    form0 = form.copy()

    nullstates = matrix_nullspace(M.E.transpose(), tol=1e-6).T  # unit vectors representing the possible nullstates

    nullspaces = []
    for nullstate in nullstates:  # reshaping the nullstates a list of (vcount x 2) arrays
        print(nullstate.shape)
        xy = nullstate.reshape((2, -1)).T
        print(xy.shape)
        nullspaces.append(xy)

    xy0 = array(form.xy())[M.free]
    delta = 1.0

    j = 0
    for nullstate in nullspaces:
        xy = xy0 + delta * nullstate

        print('Nullspace #:', j)

        for i, key in enumerate(form.vertices_where({'is_fixed': False})):
            x, y = xy[i]
            form.vertex_attribute(key, 'x', x)
            form.vertex_attribute(key, 'y', y)

        plotter = TNOPlotter(form)
        plotter.draw_form(scale_width=False)
        plotter.draw_mesh(form0, color=Color.grey())
        plotter.draw_supports()
        plotter.show()

        j += 1

# # # # ---------- EXAMPLE Barrel SAG

dx = 4
dy = 8
discretisation = [dx, dy]
xy_span = [[0, dx], [0, dy]]
tol = 1e-3

form = FormDiagram.create_ortho_form(discretisation=discretisation, xy_span=xy_span, fix='corners')
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    if abs(y - dy) < tol or y < tol:
        form.vertex_attribute(key, 'is_fixed', True)

for key in form.vertices():
    form.vertex_attribute(key, 'pz', -1.0)

for u, v in form.edges():
    if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
        form.edge_attribute((u, v), '_is_edge', False)

# add_faces = [
#     [0, 1, 3, 5, 7, 9, 11, 13, 15],
#     [44, 43, 42, 41, 40, 39, 38, 37, 36]
# ]

# for face in add_faces:
#     id = form.add_face(face)
#     form.face_attribute(id, '_is_loaded', False)

# for u, v in form.edges():
#     if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
#         form.edge_attribute((u, v), '_is_edge', False)

# delta = 0.1
# for key in form.vertices_where({'is_fixed': False}):
#     x, y, _ = form.vertex_coordinates(key)
#     dx = delta * (-1 + 2 * random_sample())
#     dy = delta * (-1 + 2 * random_sample())
#     form.vertex_attributes(key, 'xy', [x + dx, y + dy])

find_show_inds(form)

initialize_loadpath(form)
M = initialise_form(form)

# path_lp = '/Users/mricardo/compas_dev/me/inds/three_legs_lp.json'
# form.to_json(path_lp)

update_geometry(form, M)

E = M.E

find_nullspaces(form, M)

# for i, key in enumerate(form.vertices()):
#     if i in inds_column:
#         dx, dy = Et[:, i].reshape(-1, 2)

view = Viewer(form)
view.draw_thrust()
view.show()

# plotter = TNOPlotter(form)
# plotter.draw_form()
# plotter.draw_force()
# plotter.show()

# # # # # ---------- EXAMPLE 5 BARREL NO SAG

discretisation = [dx, dy]
xy_span = [[0, dx], [0, dy]]
tol = 1e-3

form = FormDiagram.create_ortho_form(discretisation=discretisation, xy_span=xy_span, fix='corners')
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    if abs(y - dy) < tol or y < tol:
        form.vertex_attribute(key, 'is_fixed', True)

for key in form.vertices():
    form.vertex_attribute(key, 'pz', -1.0)

for edge in form.edges_on_boundary():
    form.edge_attribute(edge, 'q', -10.0)

for u, v in form.edges():
    if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
        form.edge_attribute((u, v), '_is_edge', False)

# SAG ON PATTERN
equilibrium_fdm(form)

# add force diagram
# text = {key: str(key) for key in form.vertices()}

# plotter = TNOPlotter(form)
# plotter.draw_form()
# plotter.draw_vertexlabels(text=text)
# plotter.show()

# add_faces = [
#     [0, 1, 3, 5, 7, 9, 11, 13, 15],
#     [44, 43, 42, 41, 40, 39, 38, 37, 36]
# ]

# for face in add_faces:
#     id = form.add_face(face)
#     form.face_attribute(id, '_is_loaded', False)

# for u, v in form.edges():
#     if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
#         form.edge_attribute((u, v), '_is_edge', False)

# make a mess
# delta = 0.1
# for key in form.vertices_where({'is_fixed': False}):
#     x, y, _ = form.vertex_coordinates(key)
#     dx = delta * (-1 + 2 * random_sample())
#     dy = delta * (-1 + 2 * random_sample())
#     form.vertex_attributes(key, 'xy', [x + dx, y + dy])

find_show_inds(form)

initialize_loadpath(form)
M = initialise_form(form)

# path_lp = '/Users/mricardo/compas_dev/me/inds/three_legs_lp.json'
# form.to_json(path_lp)

update_geometry(form, M)

E = M.E

find_nullspaces(form, M)

update_geometry(form, M)

view = Viewer(form)
view.draw_thrust()
view.show()

# plotter = TNOPlotter(form)
# plotter.draw_form()
# plotter.draw_force()
# plotter.show()

