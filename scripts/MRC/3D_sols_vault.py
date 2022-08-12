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


discr = 14
xf = 10.0
x0 = 0.0
xc = yc = (x0 + xf)/2
xyspan = [[x0, xf], [x0, xf]]
spr_angle = 30
alpha = 1/math.cos(math.radians(spr_angle))
L = xf * alpha
Ldiff = L - xf
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=100)
cross = Shape.create_crossvault(xy_span=xyspan_shape, discretisation=discr*2)

zm = [cross.middle.vertex_attribute(key, 'z') for key in cross.middle.vertices()]
print('max, min zm', max(zm), min(zm))

# view = Viewer(form, shape=cross)
# view.draw_shape()
# view.draw_middle_shape()
# view.show()

analysis = Analysis.create_minthk_analysis(form, cross)
analysis.apply_envelope()

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x - xc > 0.1 and y - yc > 0.1:
        dXbi = normalize_vector([1, 1, 0])
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

cross = Shape.from_formdiagram_and_attributes(form)

path = '/Users/mricardo/compas_dev/me/compl_energy/crossvault/corner/cross_fd/crossvault_cross_fd_discr_14_Ecomp-linear_thk_50.0.json'
form = FormDiagram.from_json(path)

# plotter = TNOPlotter(form)
# plotter.settings['color.edges.independent'] = Color.blue()
# plotter.settings['color.vertex.supports']  = Color.red()
# plotter.settings['size.vertex'] = 6.0
# plotter.draw_form_independents()
# plotter.draw_supports()
# plotter.show()

zi = [cross.intrados.vertex_attribute(key, 'z') for key in cross.intrados.vertices()]
ze = [cross.extrados.vertex_attribute(key, 'z') for key in cross.extrados.vertices()]
zm = [cross.middle.vertex_attribute(key, 'z') for key in cross.middle.vertices()]
t = (min(ze) + min(zi))/2

print('max, min zi', max(zi), min(zi))
print('max, min ze', max(ze), min(ze))
print('max, min ze', max(zm), min(zm))

trans = Translation.from_vector(Vector(0, 0, +t))

view = Viewer(form, shape=cross)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form()
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
# view.draw_reactions()
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    base.transform(trans)
    view.draw_vector(vector=vector, base=base)
view.show()

discr = 14
xf = 10.0
x0 = 0.0
xc = yc = (x0 + xf)/2
xyspan = [[x0, xf], [x0, xf]]
spr_angle = 30
alpha = 1/math.cos(math.radians(spr_angle))
L = xf * alpha
Ldiff = L - xf
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=100)
cross = Shape.create_crossvault(xy_span=xyspan_shape, discretisation=discr*2)

analysis = Analysis.create_minthk_analysis(form, cross)
analysis.apply_envelope()

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x - xc > 0.1 and y - yc > 0.1:
        dXbi = normalize_vector([1, 1, 0])
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

cross = Shape.from_formdiagram_and_attributes(form)

path = '/Users/mricardo/compas_dev/me/compl_energy/crossvault/corner/cross_fd/crossvault_cross_fd_discr_14_Ecomp-linear_thk_50.0.json'
form = FormDiagram.from_json(path)

# plotter = TNOPlotter(form)
# plotter.settings['color.edges.independent'] = Color.blue()
# plotter.settings['color.vertex.supports']  = Color.red()
# plotter.settings['size.vertex'] = 6.0
# plotter.draw_form_independents()
# plotter.draw_supports()
# plotter.show()

zi = [cross.intrados.vertex_attribute(key, 'z') for key in cross.intrados.vertices()]
ze = [cross.extrados.vertex_attribute(key, 'z') for key in cross.extrados.vertices()]
zm = [cross.middle.vertex_attribute(key, 'z') for key in cross.middle.vertices()]
t = (min(ze) + min(zi))/2

# print('max, min zi', max(zi), min(zi))
# print('max, min ze', max(ze), min(ze))
print('aaa')
print('max, min zm', max(zm), min(zm))

# trans = Translation.from_vector(Vector(0, 0, +t))

# view = Viewer(form, shape=cross)
# view.settings['camera.show.grid'] = False
# view.settings['camera.distance'] = 35
# view.settings['camera.target'] = [5, 5, 0]
# view.settings['camera.rz'] = 45
# view.settings['camera.rx'] = 90.0
# view.settings['camera.ry'] = 0.0
# view.draw_form()
# # view.app.view.camera.rotation.y = 0.0
# view.draw_form(cull_negative=True)
# view.draw_cracks(cull_negative=True)
# # view.draw_reactions()
# view.draw_shape()
# for i in range(len(vectors_plot)):
#     vector = vectors_plot[i]
#     base = base_plot[i]
#     base.transform(trans)
#     view.draw_vector(vector=vector, base=base)
# view.show()


points = []
supports = []
edges = []

for u, v in form.edges_where({'_is_edge': True}):
    Xu = form.vertex_coordinates(u)
    Xv = form.vertex_coordinates(v)

    # if abs(Xu[0] - xc) < 1e-3 and abs(Xv[0] - xc) < 1e-3:
    if abs(Xu[1] - Xu[0]) < 1e-3 and abs(Xv[1] - Xv[0]) < 1e-3:
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

view: Viewer = Viewer(form, shape=cross)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 90.0
view.settings['camera.ry'] = 0.0
view.draw_form(edges=edges, cull_negative=True)
view.draw_cracks(points=points, cull_negative=True)
# view.draw_reactions(supports=supports, extend_reactions=False)
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()
