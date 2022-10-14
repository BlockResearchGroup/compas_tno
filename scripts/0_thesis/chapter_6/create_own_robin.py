from compas.datastructures import Mesh
from compas.geometry import distance_point_point_xy
from compas_tno import analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.problems.problems import plot_svds
from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.viewers import Viewer
from compas_tno.utilities import store_inds
from numpy import zeros
from compas_tno.problems.problems import plot_svds
from compas_tno.algorithms import constrained_smoothing

# xp, yp = 7.5, 5.0
# form = FormDiagram.create_circular_radial_form(discretisation=[4, 26])
# form.vertex_attribute(0, 'x', xp)
# form.vertex_attribute(0, 'y', yp)

# fixed = list(form.vertices_where({'is_fixed': True}))
# fixed.append(0)

# constr = {}
# for item in fixed:
#     constr[item] = form.vertex_coordinates(item)

# constrained_smoothing(form, kmax=500, constraints=constr, algorithm='area')

# plotter = TNOPlotter(form)
# plotter.draw_form(scale_width=False)
# plotter.draw_supports()
# plotter.show()

xspan = [0, 10]
yspan = [0, 10]
folder = '/Users/mricardo/compas_dev/me/pattern/singular/dome/'
prob = 'D2-mod'
mesh_file = folder + 'mesh-' + prob + '.json'
tol = None

print('File from:', mesh_file)

mesh = Mesh.from_json(mesh_file)

move_pattern_to_origin(mesh)

lines = sorted(mesh.to_lines())
meshs = Mesh.from_lines(lines)

form = FormDiagram.from_mesh(mesh)
print('form faces:', form.number_of_faces())

form.delete_boundary_edges()
form.set_boundary_supports()

# Manually find independents for this example

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
# store_inds(form)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.show()

xp, yp = 7.5, 5.0
try:
    dist = 0.01
    for key in form.vertices():
        coords = form.vertex_coordinates(key)
        print(key, coords, distance_point_point_xy(coords, [xp, yp]))
        if distance_point_point_xy(coords, [xp, yp]) < dist:
            load_node = key

    print('Loaded node is', load_node)
except:
    xo, yo = 8.0, 5.0
    dist = 0.01
    for key in form.vertices():
        coords = form.vertex_coordinates(key)
        print(key, coords, distance_point_point_xy(coords, [xp, yp]))
        if distance_point_point_xy(coords, [xo, yo]) < dist:
            load_node = key

    print('Loaded node is', load_node)

    fixed = list(form.vertices_where({'is_fixed': True}))
    fixed.append(load_node)

    constr = {}
    for item in fixed:
        constr[item] = form.vertex_coordinates(item)

    constrained_smoothing(form, kmax=500, constraints=constr, algorithm='area')

# fixed = list(form.vertices_where({'is_fixed': True}))
# fixed.append(0)

# constr = {}
# for item in fixed:
#     constr[item] = form.vertex_coordinates(item)

# constrained_smoothing(form, kmax=500, constraints=constr, algorithm='area')

# plotter = TNOPlotter(form)
# plotter.draw_form(scale_width=False)
# plotter.draw_supports()
# plotter.show()


max_load_mult = 600.0
n = form.number_of_vertices()
pzv = zeros((n, 1))
pzv[load_node] = -1.0

dome = Shape.create_dome()
dome.ro = 0.1

# analysis = Analysis.create_minthrust_analysis(form, dome,
#                                               printout=True,
#                                               plot=True,
#                                               solver='IPOPT')
analysis = Analysis.create_max_load_analysis(form, dome,
                                             load_direction=pzv,
                                             max_lambd=max_load_mult,
                                             printout=True,
                                             plot=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_optimiser_options(tol_inds=tol)
analysis.set_up_optimiser()
plot_svds(analysis.optimiser.M)
analysis.run()

view = Viewer(form)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
view.show()

# # delta = 0.2
# # for key in form.vertices_where({'is_fixed': False}):
# #     x, y, _ = form.vertex_coordinates(key)
# #     dx = delta * (-1 + 2 * random_sample())
# #     dy = delta * (-1 + 2 * random_sample())
# #     form.vertex_attributes(key, 'xy', [x + dx, y + dy])

# M = initialise_problem_general(form)

# # Find ideal tolerance

# svs = matrix_svd(M.E)
# svflip = flip(svs)

# tol = None
# maxg = 10.0
# for i in range(len(svflip) - 1):
#     dif = svflip[i+1] - svflip[i]
#     pc = dif/svflip[i]
#     print(i, pc, svflip[i+1], svflip[i])
#     if pc > maxg:
#         tol = (svflip[i+1] + svflip[i])/2
#         break

# print('Suggested tolerance:', tol)

# start_time = time.time()

# # tol = None
# ind = find_independents_forward(M.E, tol=tol)
# # ind = find_independents_backward(M.E, tol=tol)

# elapsed_time = time.time() - start_time
# print('Elapsed Time:', elapsed_time)

# index_uv = form.index_uv()
# ind_edges = [index_uv[i] for i in ind]

# store_inds(form, ind_edges=ind_edges)

# #     print('Max / min svds', i, max(svds), min(svds))

# adapt_problem_to_fixed_diagram(M, form, printout=True)

# # plotter = TNOPlotter(form)
# # plotter.draw_form_independents()
# # plotter.draw_supports()
# # plotter.show()

# plot_svds(M)

# # print(matrix_svd(M.E[:, M.dep]))

# # plotter = TNOPlotter(form)
# # plotter.draw_form_independents()
# # plotter.formartist.draw_vertexlabels()  # draw_vertices(facetext={key: key for key in form.vertices()})
# # plotter.draw_supports()
# # plotter.show()

# dome = Shape.create_dome()

# analysis = Analysis.create_minthrust_analysis(form, dome, printout=True)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.set_up_optimiser()
# analysis.run()

# plotter = TNOPlotter(form, form_base=form_)
# plotter.draw_base_form()
# plotter.show_solution()
