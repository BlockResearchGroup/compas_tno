import math
import compas_tno

from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import move_pattern_to_origin

from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Scale
from compas.datastructures import Mesh

from numpy import array

import os

xf = 10.0
x0 = 0.0
yf = 10.0
y0 = 0.0

spr_angle = 30.0
discretisation = 14

displ_type = 'two'

for thk in [0.50]:

    pavillion = Shape.create_pavillionvault(thk=thk, spr_angle=spr_angle, expanded=True)
    pavillion.ro = 20.0

    # form = FormDiagram.create_ortho_form(fix='all', discretisation=discretisation)
    form = FormDiagram.create_cross_form(fix='all', discretisation=discretisation)
    # form = FormDiagram.create_cross_with_diagonal(fix='all', discretisation=discretisation)

    coef = 1/math.cos(math.radians(spr_angle))

    for key in form.vertices_where({'is_fixed': True}):
        x, y, _ = form.vertex_coordinates(key)
        bx = 0
        by = 0
        if abs(x - xf) < 1e-3:
            bx = + coef * thk / 2
        if abs(x - x0) < 1e-3:
            bx = - coef * thk / 2
        if abs(y - yf) < 1e-3:
            by = + coef * thk / 2
        if abs(y - y0) < 1e-3:
            by = - coef * thk / 2
        form.vertex_attribute(key, 'b', [bx, by])

    # path = '/Users/mricardo/compas_dev/me/pattern/equidistant/equidistant_discr_14.json'
    # mesh = Mesh.from_json(path)

    # move_pattern_to_origin(mesh, corners=[[0.001, 0.001], [9.999, 0.001], [9.999, 9.999], [0.001, 9.999]])

    # form = FormDiagram.from_mesh(mesh)
    # form.parameters['type'] = 'equidistant'

    # form.delete_boundary_edges()
    # form.set_boundary_supports()

    # pl = TNOPlotter(form)
    # pl.draw_form(scale_width=False)
    # pl.draw_supports()
    # pl.show()

    # view = Viewer(form, shape=pavillion)
    # view.draw_thrust()
    # view.draw_shape()
    # view.draw_b_constraint()
    # view.show()

    vector_supports = []
    vectors_plot = []
    base_plot = []

    for key in form.vertices_where({'is_fixed': True}):
        x, y, z = form.vertex_coordinates(key)
        dXbi = [0, 0, 0]
        if abs(x - xf) < 0.01:
            dXbi = [1, 0, 0]
            vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
            base_plot.append(Point(x, y, z - 0.2))
        # if abs(x - x0) < 0.01:
        #     dXbi = [-1, 0, 0]
        #     vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
        #     base_plot.append(Point(x, y, z - 0.2))

        vector_supports.append(dXbi)

    dXb = array(vector_supports)

    constraints = ['funicular', 'envelope', 'reac_bounds']
    # features = ['fixed', 'sym']
    features = ['fixed']
    axis_sym = [[[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]], [[5.0, 0.0, 0.0], [5.0, 10.0, 0.0]]]
    starting = 'loadpath'
    # starting = 'current'
    # constraints = ['funicular', 'envelope']

    # form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/form.json')

    analysis = Analysis.create_compl_energy_analysis(form,
                                                     pavillion,
                                                     printout=True,
                                                     plot=True,
                                                     support_displacement=dXb,
                                                     max_iter=1000,
                                                     solver='IPOPT',
                                                     starting_point=starting)

    analysis.optimiser.set_constraints(constraints)
    analysis.optimiser.set_features(features)
    analysis.optimiser.set_axis_symmetry(axis_sym)
    analysis.apply_selfweight()
    analysis.apply_envelope()

    # for key in form.vertices_where({'is_fixed': True}):
    #     lb = form.vertex_attribute(key, 'lb')
    #     form.vertex_attribute(key, 'lb', lb - 0.5)

    # analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()

    from compas_tno.problems.problems import plot_svds
    M = analysis.optimiser.M
    plot_svds(M)

    analysis.run()

    # address = compas_tno.get('form.json')
    # form.to_json(address)
    # print('Saved Last iteration to:', address)

    # folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open', form.parameters['type']) + '/'
    # os.makedirs(folder, exist_ok=True)
    # title = pavillion.datashape['type'] + '_' + form.parameters['type'] + '_discr_' + str(discretisation)
    # if spr_angle:
    #     title = pavillion.datashape['type'] + '_' + form.parameters['type'] + '_discr_' + str(discretisation) + '_spr_' + str(spr_angle)
    # address = folder + title + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
    # if analysis.optimiser.exitflag == 0:
    #     print('Saving form to:', address)
    #     form.to_json(address)

plotter = TNOPlotter(form)
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    plotter.draw_vector(vector=vector, base=base)
plotter.draw_form()
plotter.draw_supports()
plotter.draw_cracks()
plotter.show()

view = Viewer(form, shape=pavillion)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()

view = Viewer(form, shape=pavillion)
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()
