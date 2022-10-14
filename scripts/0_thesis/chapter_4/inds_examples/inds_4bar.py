from compas_plotters import Plotter
from compas.geometry import Line
from compas.datastructures import Mesh
from compas.geometry import distance_point_point_xy
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.problems import initialise_problem_general
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.problems.problems import plot_svds
from numpy.random import random_sample

pts = [[0, 0], [-2, 0], [0, 2], [2, 0], [0, -2]]
fixed_pts = [[-2, 0], [0, 2], [2, 0], [0, -2]]

lines = [Line(pts[0], pts[1]),
         Line(pts[0], pts[2]),
         Line(pts[0], pts[3]),
         Line(pts[0], pts[4]),
         ]

mesh = Mesh.from_lines(lines)
form = FormDiagram.from_mesh(mesh)

for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    for sup in fixed_pts:
        if distance_point_point_xy([x, y], sup) < 1e-3:
            form.vertex_attribute(key, 'is_fixed', True)

# UNcomment to check the case with random-ness
# delta = 0.5
# for key in form.vertices_where({'is_fixed': False}):
#     x, y, _ = form.vertex_coordinates(key)
#     dx = delta * (-1 + 2 * random_sample())
#     dy = delta * (-1 + 2 * random_sample())
#     form.vertex_attributes(key, 'xy', [x + dx, y + dy])

M = initialise_problem_general(form)
adapt_problem_to_fixed_diagram(M, form, printout=True)

plot_svds(M)

plotter = TNOPlotter(form)
plotter.draw_form_independents()
plotter.draw_supports()
plotter.show()
