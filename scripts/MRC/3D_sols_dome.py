from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import apply_bounds_reactions
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape

from compas.geometry import Vector
from compas.geometry import Point
from compas.colors import Color

from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import normalize_vector
from compas.geometry import dot_vectors

from numpy import array, linspace
from compas_tno.shapes import Shape, rectangular_topology
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.shapes.pavillionvault import pavillionvault_middle_update, pavillionvault_ub_lb_update

from numpy import array
import math

import os


x0, xf = 0.0, 10.0
xc = (x0 + xf)/2
yc = xc

thk = 0.50
form = FormDiagram.create_circular_radial_form()

vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)

    # Outward displacement
    # dXbi = normalize_vector([(x - xc), (y - yc), 0.0])
    # vectors_plot.append(Vector(*dXbi))
    # base_plot.append(Point(x, y, z))

    # Split Displacement
    dXbi = [0, 0, 0]
    if x - xc > 0.1:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))
    if x - xc < -0.1:
        dXbi = [-1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

dome = Shape.create_dome_polar(discretisation=[50, 50])

# Traditional form diagram
paths = ['/Users/mricardo/compas_dev/me/compl_energy/dome/outwards/radial_fd/dome_radial_fd_discr_[16, 20]_Ecomp-linear_thk_50.0.json']

# Non traditional form diagram
# paths = ['/Users/mricardo/compas_dev/me/compl_energy/dome/split/radial_fd/dome_radial_fd_discr_[16, 20]_Ecomp-linear_thk_50.0.json']
paths = []
for ni in [3]:  # [4, 2], [3, 6], [2, 5], [1, 4]
    for ns in [6]:
        pt = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/parametric/form_sigularity_ni_{0}_ns_{1}.json/circular_sigularity_ni_{0}_ns_{1}.json_Ecomp-linear_thk_50.0.json'.format(ni, ns)
        paths.append(pt)

for path in paths:
    form = FormDiagram.from_json(path)

    i = 0
    dic = {}
    for key in form.vertices_where({'is_fixed': True}):
        rx, ry, rz = form.vertex_attributes(key, ['_rx', '_ry', '_rz'])
        print(i, key, rx, ry, rz)
        dic[key] = i
        i += 1

    plotter = TNOPlotter(form)
    plotter.draw_form()
    plotter.draw_cracks()
    # plotter.draw_form_independents()
    plotter.draw_vertexlabels(dic, 20)
    plotter.draw_supports()
    plotter.show()

    plotter = TNOPlotter(form)
    plotter.settings['color.edges.independent'] = Color.blue()
    plotter.settings['color.vertex.supports']  = Color.red()
    plotter.settings['size.vertex'] = 6.0
    plotter.draw_form_independents()
    plotter.draw_supports()
    plotter.show()


    view = Viewer(form, shape=dome)
    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['camera.rz'] = 45
    view.settings['camera.rx'] = 60
    view.draw_form(cull_negative=True)
    view.draw_cracks(cull_negative=True)
    view.draw_reactions(extend_reactions=True)
    view.draw_shape()
    for i in range(len(vectors_plot)):
        vector = vectors_plot[i]
        base = base_plot[i]
        view.draw_vector(vector=vector, base=base)
    view.show()

    points = []
    supports = []
    edges = []

    for u, v in form.edges_where({'_is_edge': True}):
        Xu = form.vertex_coordinates(u)
        Xv = form.vertex_coordinates(v)

        # if abs(Xu[0] - xc) < 1e-3 and abs(Xv[0] - xc) < 1e-3:
        if abs(Xu[1] - yc) < 1e-3 and abs(Xv[1] - yc) < 1e-3:
            edges.append((u, v))
            if u not in points:
                points.append(u)
            if v not in points:
                points.append(v)
            if form.vertex_attributes(u, 'is_fixed'):
                supports.append(u)
            if form.vertex_attributes(v, 'is_fixed'):
                supports.append(v)

    # delete_faces = []
    # for face in pavillion.intrados.faces():
    #     Xf = pavillion.intrados.face_centroid(face)
    #     if Xf[1] < yc:
    #         delete_faces.append(face)

    # for face in delete_faces:
    #     pavillion.intrados.delete_face(face)
    #     pavillion.extrados.delete_face(face)

    view: Viewer = Viewer(form, shape=dome)
    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 15
    view.settings['size.vertex'] = 16.0
    view.settings['size.edge.max_thickness'] = 12.0
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['camera.rz'] = 0
    view.settings['camera.rx'] = 90
    view.draw_form(edges=edges, cull_negative=True)
    view.draw_cracks(points=points, cull_negative=True)
    view.draw_reactions(supports=supports, extend_reactions=True)
    view.draw_shape()
    for i in range(len(vectors_plot)):
        vector = vectors_plot[i]
        base = base_plot[i]
        view.draw_vector(vector=vector, base=base)
    view.show()

    view = Viewer(form, shape=dome)
    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['camera.rz'] = 45
    view.settings['camera.rx'] = 60
    # view.draw_form(cull_negative=True)
    # view.draw_cracks(cull_negative=True)
    # view.draw_reactions(extend_reactions=True)
    view.draw_shape()
    # for i in range(len(vectors_plot)):
    #     vector = vectors_plot[i]
    #     base = base_plot[i]
    #     view.draw_vector(vector=vector, base=base)
    view.show()

    plotter = TNOPlotter(form, dome, figsize=(16,5))
    plotter.draw_form()
    plotter.draw_cracks()
    plotter.draw_supports()
    plotter.draw_force()
    plotter.show()

    force = plotter.force
    force_path = path.split('.')[0] + '_force.json'
    force.to_json(force_path)
    print('Force Saved:', force_path)

