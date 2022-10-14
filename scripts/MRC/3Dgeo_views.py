from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape, rectangular_topology
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.shapes.pavillionvault import pavillionvault_middle_update, pavillionvault_ub_lb_update
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import apply_bounds_reactions
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape

from compas.geometry import Vector
from compas.geometry import Point

from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import normalize_vector

from numpy import array, linspace
from compas_tno.shapes import Shape, rectangular_topology
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.shapes.pavillionvault import pavillionvault_middle_update, pavillionvault_ub_lb_update
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
thk = 0.50
discretisation = 14
xc = yc = (xspan[1] - xspan[0])/2
x0 = y0 = xspan[0]
x1 = y1 = xspan[1]
form = FormDiagram.create_cross_with_diagonal(fix='all', discretisation=discretisation)
xy_span = [xspan, yspan]

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if abs(x - xspan[1]) < 1e-3:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

spr_angle = 30.0

thk = 0.50

pavillion = pavillion_nice(thk, xy_span, spr_angle)

view = Viewer(form, shape=pavillion)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()

view = Viewer(form, shape=pavillion)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_shape()
    # for i in rang{{{e(len(vectors_plot)):
    #     vector = vectors_plot[i]
    #     base = base_plot[i]
    #     view.dr}}}aw_vector(vector=vector, base=base)
view.show()


# OLD PAVILLION

# pavillion_view = Shape.create_pavillionvault(thk=thk, t=0.0, expanded=True)

# intra = pavillion_view.intrados
# extra = pavillion_view.extrados
# middle = pavillion_view.middle

# scale = Scale.from_factors([1.0, 1.0, 3.5/5.0])  # change apex height to 3.5
# intra.transform(scale)
# extra.transform(scale)
# middle.transform(scale)

# pavillion_nice = Shape.from_meshes(intra, extra, middle)

# view = Viewer(form, shape=pavillion_nice)
# view.settings['camera.show.grid'] = False
# view.settings['camera.distance'] = 35
# view.settings['camera.target'] = [5, 5, 0]
# view.settings['camera.rz'] = 45
# view.settings['camera.rx'] = 60
# view.draw_shape()
# for i in range(len(vectors_plot)):
#     vector = vectors_plot[i]
#     base = base_plot[i]
#     view.draw_vector(vector=vector, base=base)
# view.show()


#------------ DOME PLOT

form = FormDiagram.create_circular_radial_form()

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x - xc > 0.1:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))
    if x - xc < -0.1:
        dXbi = [-1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

dome = Shape.create_dome_polar(discretisation=[50, 50])

view = Viewer(form, shape=dome)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()


# ---------- CROSS

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

zi = [cross.intrados.vertex_attribute(key, 'z') for key in cross.intrados.vertices()]
ze = [cross.extrados.vertex_attribute(key, 'z') for key in cross.extrados.vertices()]
t = (min(ze) + min(zi))/2

print('max, min zi', max(zi), min(zi))
print('max, min ze', max(ze), min(ze))

trans = Translation.from_vector(Vector(0, 0, -t))
cross.intrados.transform(trans)
cross.extrados.transform(trans)

view = Viewer(form, shape=cross)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()

view = Viewer(form, shape=cross)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_shape()
# for i in range(len(vectors_plot)):
#     vector = vectors_plot[i]
#     base = base_plot[i]
#     view.draw_vector(vector=vector, base=base)
view.show()
