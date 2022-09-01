from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities.form import move_pattern_to_origin
from compas_tno.viewers import Viewer
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.problems.problems import plot_svds

from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Line
from compas.geometry import norm_vector
from compas.geometry import normalize_vector
from compas.geometry import distance_point_point_xy
from compas.colors import Color

import compas_tno

from numpy import array

import math

discr = 14
xf = yf = 10.0
x0 = y0 = 0.0
xc = yc = (x0 + xf)/2
xyspan = [[x0, xf], [x0, xf]]
spr_angle = 30
alpha = 1/math.cos(math.radians(spr_angle))
L = xf * alpha
Ldiff = L - xf
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

path = '/Users/mricardo/compas_dev/me/pattern/crossvault/form-cross+crack.json'

form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=discr)
# form = FormDiagram.create_fan_form(xy_span=xyspan, discretisation=discr)

# form = FormDiagram.from_json(path)
# move_pattern_to_origin(form, fix_corners=True)

plotter: TNOPlotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.show()

shape = Shape.create_crossvault(xy_span=xyspan_shape, discretisation=discr*2)

vector_supports = []
lines = []
plots_vectors = []
for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x > xc and y < yc:
        # dXbi = normalize_vector([1, -1, 0])  # outwards
        dXbi = normalize_vector([-1, +1, 0])  # inwards
    if norm_vector(dXbi) > 1e-3:
        start = [x, y, z]
        end = [x + dXbi[0], y + dXbi[1], z + dXbi[2]]
        lines.append(Line(start, end))
        plots_vectors.append([Point(*start), Vector(*dXbi)])
    vector_supports.append(dXbi)

dXb = array(vector_supports)

plotter: TNOPlotter = TNOPlotter(form, shape=shape)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
for base, vector in plots_vectors:
    plotter.draw_vector(vector, base)
plotter.show()

analysis: Analysis = Analysis.create_compl_energy_analysis(form, shape, printout=True, support_displacement=dXb)
# analysis.optimiser.settings['tol_inds'] = 1e-6
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()

plot_svds(analysis.optimiser.M)

plotter: TNOPlotter = TNOPlotter(form, shape=shape)
plotter.draw_form_independents()
plotter.draw_supports()
for base, vector in plots_vectors:
    plotter.draw_vector(vector, base)
plotter.show()

analysis.run()

path = compas_tno.get('form.json')
form.to_json(path)
print('Save Form to:', path)

print('SWT:', form.lumped_swt())

plotter: TNOPlotter = TNOPlotter(form, shape=shape)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.draw_cracks()
for base, vector in plots_vectors:
    plotter.draw_vector(vector, base)
plotter.show()

view: Viewer = Viewer(form, show_grid=False)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
# view.draw_reactions()
for base, vector in plots_vectors:
    print(base, vector)
    view.draw_vector(vector, base, color=Color.black())
view.show()
