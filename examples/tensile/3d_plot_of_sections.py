from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas.geometry import Vector
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Polygon
from compas.geometry import Rotation
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape

delta = 1.0
span = 10.0
xspan = yspan = [0.0, span]
x0 = y0 = 0.0
x1 = y1 = span
tol = 1e-4
xspan_vault = yspan_vault = [- delta, span + delta]
thk = 0.5
discretisation = 10

data_vault = {
    'type': 'crossvault',
    'xy_span': [xspan_vault, yspan_vault],
    'thk': thk,
    'discretisation': discretisation,
    't': 0.0
}

data_diagram = {
    'type': 'cross_fd',
    'xy_span': [xspan, yspan],
    'thk': thk,
    'discretisation': discretisation,
    'fix': 'corners'
}

vault = Shape.from_library(data_vault)
form = FormDiagram.from_library(data_diagram)

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)
vault = Shape.from_formdiagram_and_attributes(form)

form.vertices_attribute('wx', span/discretisation)
form.vertices_attribute('wy', span/discretisation)

viewer1 = Viewer(form, vault)
viewer1.settings['camera.show.grid'] = False
viewer1.settings['camera.distance'] = 35
viewer1.draw_shape()

viewer2 = Viewer(form, vault)
viewer2.settings['camera.show.grid'] = False
viewer2.settings['camera.distance'] = 35
viewer2.draw_shape()

viewer3 = Viewer(form, vault)
viewer3.settings['camera.show.grid'] = False
viewer3.settings['camera.distance'] = 35
viewer3.draw_shape()

i = Vector(1.0, 0, 0)
j = Vector(0, 1.0, 0)
k = Vector(0, 0, 1.0)

for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    z = form.vertex_attribute(key, 'target')
    ub, lb = form.vertex_attribute(key, 'ub'), form.vertex_attribute(key, 'lb')
    t = ub - lb
    wxi = form.vertex_attribute(key, 'wx')
    wyi = form.vertex_attribute(key, 'wy')
    line1 = Line([x - wxi/2, y, z], [x + wxi/2, y, z])
    line2 = Line([x, y - wyi/2, z], [x, y + wyi/2, z])
    normal = vault.middle.vertex_normal(key)

    # compute "diaedric angles"
    angle_x = i.angle(normal)
    angle_y = j.angle(normal)
    angle_z = k.angle(normal)

    # create a polygon centered in the middle surface
    polygon1 = Polygon([
        line1.start + Point(0, 0, -t/2),
        line1.end + Point(0, 0, -t/2),
        line1.end + Point(0, 0, +t/2),
        line1.start + Point(0, 0, +t/2)
    ])

    polygon2 = Polygon([
        line2.start + Point(0, 0, -t/2),
        line2.end + Point(0, 0, -t/2),
        line2.end + Point(0, 0, +t/2),
        line2.start + Point(0, 0, +t/2)
        ])

    if y <= y0 + (y1 - y0)/(x1 - x0) * (x - x0) + tol and y >= y1 - (y1 - y0)/(x1 - x0) * (x - x0) - tol:  # Q1
        if y < (y1 + y0)/2:
            vector_rotation = [+angle_z, 0, 0]
        else:
            vector_rotation = [-angle_z, 0, 0]
    elif y >= y0 + (y1 - y0)/(x1 - x0) * (x - x0) - tol and y >= y1 - (y1 - y0)/(x1 - x0) * (x - x0) - tol:  # Q3
        if x < (x1 + x0)/2:
            vector_rotation = [0, -angle_z, 0]
        else:
            vector_rotation = [0, +angle_z, 0]
    elif y >= y0 + (y1 - y0)/(x1 - x0) * (x - x0) - tol and y <= y1 - (y1 - y0)/(x1 - x0) * (x - x0) + tol:  # Q2
        if y < (y1 + y0)/2:
            vector_rotation = [+angle_z, 0, 0]
        else:
            vector_rotation = [-angle_z, 0, 0]
    elif y <= y0 + (y1 - y0)/(x1 - x0) * (x - x0) + tol and y <= y1 - (y1 - y0)/(x1 - x0) * (x - x0) + tol:  # Q4
        if x < (x1 + x0)/2:
            vector_rotation = [0, -angle_z, 0]
        else:
            vector_rotation = [0, +angle_z, 0]
    else:
        pass

    transf = Rotation.from_axis_angle_vector(axis_angle_vector=vector_rotation, point=[x, y, z])
    polygon1 = polygon1.transformed(transf)
    viewer1.app.add(polygon1, linewidth=5.0)

    transf = Rotation.from_axis_angle_vector(axis_angle_vector=vector_rotation, point=[x, y, z])
    polygon2 = polygon2.transformed(transf)
    viewer2.app.add(polygon2, linewidth=5.0)

    viewer3.app.add(polygon1, linewidth=5.0)
    viewer3.app.add(polygon2, linewidth=5.0)

viewer1.show()
viewer2.show()
viewer3.show()
