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
from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.utilities import store_inds
from compas_tno.algorithms.independents import count_inds, find_independents_backward, find_independents_forward
from compas_tno.algorithms.independents import matrix_svd

from numpy.linalg import matrix_rank
from numpy import flip

import time

# def find_tolerance(M):

#     s = matrix_svd(M.E)
#     if min(s) < 1e-6

xspan = [0, 10]
yspan = [0, 10]
folder = '/Users/mricardo/compas_dev/me/pattern/singular/dome/'
prob = 'A2-sym'
mesh_file = folder + 'mesh-' + prob + '.json'

mesh = Mesh.from_json(mesh_file)

move_pattern_to_origin(mesh)

form = FormDiagram.from_mesh(mesh)
print('form faces:', form.number_of_faces())

form.delete_boundary_edges()
form.set_boundary_supports()

form_ = form.copy()

# # Manually find independents for this example

# fixed = list(form.vertices_where({'is_fixed': True}))
# edges = list(form.edges_where({'_is_edge': True}))
# # inds manually
# i = 0
# for u, v in edges:
#     if i < len(fixed) - 2:
#         if u in fixed or v in fixed:
#             form.edge_attribute((u, v), 'is_ind', True)
#             i = i + 1

# edges_hoops = [(78, 75), (49, 79), (9, 141), (99, 152), (138, 56)]
# for edge in edges_hoops:
#     if not edge in edges:
#         edge = (edge[1], edge[0])
#     form.edge_attribute(edge, 'is_ind', True)

# from compas_tno.utilities import store_inds
# store_inds(form)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.show()

# delta = 0.2
# for key in form.vertices_where({'is_fixed': False}):
#     x, y, _ = form.vertex_coordinates(key)
#     dx = delta * (-1 + 2 * random_sample())
#     dy = delta * (-1 + 2 * random_sample())
#     form.vertex_attributes(key, 'xy', [x + dx, y + dy])

M = initialise_problem_general(form)

# Find ideal tolerance

svs = matrix_svd(M.E)
svflip = flip(svs)

tol = None
maxg = 10.0
for i in range(len(svflip) - 1):
    dif = svflip[i+1] - svflip[i]
    pc = dif/svflip[i]
    print(i, pc, svflip[i+1], svflip[i])
    if pc > maxg:
        tol = (svflip[i+1] + svflip[i])/2
        break

print('Suggested tolerance:', tol)

start_time = time.time()

# tol = None
ind = find_independents_forward(M.E, tol=tol)
# ind = find_independents_backward(M.E, tol=tol)

elapsed_time = time.time() - start_time
print('Elapsed Time:', elapsed_time)

index_uv = form.index_uv()
ind_edges = [index_uv[i] for i in ind]

store_inds(form, ind_edges=ind_edges)

#     print('Max / min svds', i, max(svds), min(svds))

adapt_problem_to_fixed_diagram(M, form, printout=True)

# plotter = TNOPlotter(form)
# plotter.draw_form_independents()
# plotter.draw_supports()
# plotter.show()

plot_svds(M)

# print(matrix_svd(M.E[:, M.dep]))

# plotter = TNOPlotter(form)
# plotter.draw_form_independents()
# plotter.formartist.draw_vertexlabels()  # draw_vertices(facetext={key: key for key in form.vertices()})
# plotter.draw_supports()
# plotter.show()

dome = Shape.create_dome()

analysis = Analysis.create_minthrust_analysis(form, dome, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

plotter = TNOPlotter(form, form_base=form_)
plotter.draw_base_form()
plotter.show_solution()
