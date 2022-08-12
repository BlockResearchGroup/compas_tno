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

r = 1.0
diagonal = False
if diagonal:
    form = FormDiagram.create_circular_radial_form(discretisation=[2, 4], radius=r, diagonal=True, partial_diagonal='right')
else:
    form = FormDiagram.create_circular_radial_form(discretisation=[2, 4], radius=r)

total = 0
total_length = 0.0
for edge in form.edges_where({'_is_edge': True}):
    el = form.edge_length(*edge)
    total_length += el
    total += 1
average = total_length/total
print('Edge average length:', average)

deltapc = 0.00000000000001

delta = average * deltapc
# delta = 0.001
# delta = 0.0001  # <- works
print('Delta/Length:', delta/average)

for key in form.vertices_where({'is_fixed': False}):
    x, y, _ = form.vertex_coordinates(key)
    dx = delta * (-1 + 2 * random_sample())
    dy = delta * (-1 + 2 * random_sample())
    form.vertex_attributes(key, 'xy', [x + dx, y + dy])

M = initialise_problem_general(form)
adapt_problem_to_fixed_diagram(M, form, printout=True)

plot_svds(M)
n, m = M.E.shape

from scipy.optimize import nnls
from numpy import zeros

sols, norm = nnls(M.E, zeros((n)))
print(sols)
print(norm)

plotter = TNOPlotter(form)
plotter.draw_form_independents()
plotter.draw_supports()
plotter.zoom_extents()

plotter.show()
