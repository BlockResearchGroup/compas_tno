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

from numpy import array, linspace
from compas_tno.shapes import Shape, rectangular_topology
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.shapes.pavillionvault import pavillionvault_middle_update, pavillionvault_ub_lb_update

from numpy import array
import math

import os


def pavillion_nice(thk, xy_span, spr_angle):

    [[x0, x1],[y0, y1]] = xy_span

    discretisation = 50

    x = linspace(x0, x1, num=discretisation+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0, y1, num=discretisation+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)

    xi, yi, faces_i = rectangular_topology(x, y)

    zt = pavillionvault_middle_update(xi, yi, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)
    xyzt = array([xi, yi, zt.flatten()]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    x = linspace(x0 - thk/2 / math.cos(math.radians(spr_angle)), x1 + thk/2 / math.cos(math.radians(spr_angle)), num=discretisation+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0 - thk/2 / math.cos(math.radians(spr_angle)), y1 + thk/2 / math.cos(math.radians(spr_angle)), num=discretisation+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)
    xi, yi, faces_i = rectangular_topology(x, y)
    zub, _ = pavillionvault_ub_lb_update(xi, yi, thk, 0.0, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)
    xyzub = array([xi, yi, zub.flatten()]).transpose()
    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)

    x = linspace(x0 + thk/2 / math.cos(math.radians(spr_angle)), x1 - thk/2 / math.cos(math.radians(spr_angle)), num=discretisation+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0 + thk/2 / math.cos(math.radians(spr_angle)), y1 - thk/2 / math.cos(math.radians(spr_angle)), num=discretisation+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)
    xi, yi, faces_i = rectangular_topology(x, y)
    _, zlb = pavillionvault_ub_lb_update(xi, yi, thk, 0.0, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)
    xyzlb = array([xi, yi, zlb.flatten()]).transpose()
    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)

    pavillion = Shape.from_meshes(intrados, extrados, middle)

    return pavillion


xspan = yspan = [0, 10.0]
discretisation = 14
xc = yc = (xspan[1] - xspan[0])/2
xy_span = [xspan, yspan]

# ------------ PAVILLION PLOT

spr_angle = 30.0

form_type = 'cross_fd'
shape_type = 'pavillionvault'
objective = 'Ecomp-linear'

for thk in [0.50]:

    folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open', form_type) + '/'
    title = shape_type + '_' + form_type + '_discr_' + str(discretisation) + '_spr_' + str(spr_angle)
    address = folder + title + '_' + objective + '_thk_' + str(100*thk) + '.json'

    form = FormDiagram.from_json(address)
    print('Loaded form diagram from:', address)

    vector_supports = []
    vectors_plot = []
    base_plot = []

    for key in form.vertices_where({'is_fixed': True}):
        x, y, z = form.vertex_coordinates(key)
        dXbi = [0, 0, 0]
        if abs(x - xspan[1]) < 1e-3:
            dXbi = [1, 0, 0]
            vectors_plot.append(Vector(*dXbi))
            base_plot.append(Point(x, y, 0.0))

        vector_supports.append(dXbi)

    pavillion = pavillion_nice(thk, xy_span, spr_angle)

    zmax = 0.0
    for key in pavillion.middle.vertices():
        _, _, z = pavillion.middle.vertex_coordinates(key)
        if z > zmax:
            zmax = z
    print('Maximum Z:', zmax)

    # plot: TNOPlotter = TNOPlotter(form)
    # plot.draw_form(scale_width=False, color=Color.black())
    # plot.draw_supports(color=Color.red())
    # plot.show()

    # view: Viewer = Viewer(form, shape=pavillion)
    # view.settings['camera.show.grid'] = False
    # view.settings['camera.distance'] = 35
    # view.settings['size.vertex'] = 16.0
    # view.settings['size.edge.max_thickness'] = 12.0
    # view.settings['camera.target'] = [5, 5, 0]
    # view.settings['camera.rz'] = 45
    # view.settings['camera.rx'] = 60
    # view.draw_form(cull_negative=True)
    # view.draw_cracks(cull_negative=True)
    # view.draw_reactions(extend_reactions=True)
    # view.draw_shape()
    # for i in range(len(vectors_plot)):
    #     vector = vectors_plot[i]
    #     base = base_plot[i]
    #     view.draw_vector(vector=vector, base=base)
    # view.show()

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

    delete_faces = []
    for face in pavillion.intrados.faces():
        Xf = pavillion.intrados.face_centroid(face)
        if Xf[1] < yc:
            delete_faces.append(face)

    for face in delete_faces:
        pavillion.intrados.delete_face(face)
        pavillion.extrados.delete_face(face)

    view: Viewer = Viewer(form, shape=pavillion)
    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
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

    # plot: TNOPlotter = TNOPlotter(form, shape=pavillion)
    # plot.draw_form_xz(edges=edges, scale_width=True)
    # plot.draw_shape_xz()
    # plot.draw_cracks(points=points)
    # # plot.draw_supports(color=Color.red())
    # plot.show()


