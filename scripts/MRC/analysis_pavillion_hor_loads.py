import math
import compas_tno

from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.utilities import slide_diagram

from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Scale
from compas.datastructures import Mesh

from compas_tno.utilities import apply_horizontal_multiplier

from numpy import array

import os

xf = 10.0
x0 = 0.0
yf = 10.0
y0 = 0.0

xc = (xf - x0)/2
yc = (yf - y0)/2

spr_angle = 30.0
discretisation = 14


displ_type = 'two'

lambd0 = 0.3

slide = True
delta = 0.1

for thk in [0.50]:

    pavillion = Shape.create_pavillionvault(thk=thk, spr_angle=spr_angle, expanded=True, t= 0.1)

    # form = FormDiagram.create_ortho_form(fix='all', discretisation=discretisation)
    form = FormDiagram.create_cross_form(fix='all', discretisation=discretisation)
    # form = FormDiagram.create_cross_with_diagonal(fix='all', discretisation=discretisation)

    if slide:
        slide_diagram(form, delta=-delta)

    # delta = 10/14/2

    # yc = xc = (yf - y0)/2
    # for vertex in form.vertices_where({'is_fixed': False}):
    #     x, y, _ = form.vertex_coordinates(vertex)
    #     dy = min(y - y0, yf - y)
    #     if abs(dy) > 1e-3:
    #         dx = delta * (1 - ((dy - yc)/yc)**2)
    #         if (x - xc) < - 0.01:
    #             form.vertex_attribute(vertex, 'x', x + dx)
    #         elif (x - xc) > 0.01:
    #             form.vertex_attribute(vertex, 'x', x - dx)

    # plotter = TNOPlotter(form)
    # plotter.draw_form(scale_width=False)
    # plotter.draw_supports()
    # plotter.show()

    coef = 1/math.cos(math.radians(spr_angle))

    for key in form.vertices_where({'is_fixed': True}):
        x, y, _ = form.vertex_coordinates(key)
        bx = thk
        by = thk
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
            dXbi = [+1, 0, 0]
            vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
            base_plot.append(Point(x, y, z - 0.2))

        # if abs(x - xf) < 0.01 and y - yc > 0.0:
        #     dXbi = [+1, 0, 0]
        #     vectors_plot.append(Vector(dXbi[0], dXbi[1], 0.0))
        #     base_plot.append(Point(x, y, z - 0.2))
        # if abs(y - yf) < 0.01 and x - xc > 0.0:
        #     dXbi = [0, +1, 0]
        #     vectors_plot.append(Vector(dXbi[0], dXbi[1], 0.0))
        #     base_plot.append(Point(x, y, z - 0.2))

        # if abs(x - x0) < 0.01:
        #     dXbi = [-1, 0, 0]
        #     vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
        #     base_plot.append(Point(x, y, z - 0.2))

        vector_supports.append(dXbi)

    dXb = array(vector_supports)

    plotter: TNOPlotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_supports()
    for i in range(len(vectors_plot)):
        vector = vectors_plot[i]
        base = base_plot[i]
        plotter.draw_vector(vector=vector, base=base)
    plotter.show()

    constraints = ['funicular', 'envelope', 'reac_bounds']
    features = ['fixed', 'sym']
    # features = ['fixed']
    # axis_sym = [[[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]], [[5.0, 0.0, 0.0], [5.0, 10.0, 0.0]]]
    axis_sym = [[[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]]]
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
    # analysis.optimiser.set_features(features)
    # analysis.optimiser.set_axis_symmetry(axis_sym)
    analysis.apply_selfweight()
    apply_horizontal_multiplier(form, lambd=lambd0)
    analysis.apply_envelope()

    print('SWT:', form.lumped_swt())

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

    folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open', form.parameters['type'], 'hor_loads') + '/'
    os.makedirs(folder, exist_ok=True)
    title = 'inverse_displ_' + pavillion.datashape['type'] + '_' + form.parameters['type'] + '_discr_' + str(discretisation) + '_lambdh_' + str(lambd0) + '_'
    if spr_angle:
        title = title + '_spr_' + str(spr_angle)
    if slide:
        title = title + 'sliding_diagram_' + str(delta)
    address = folder + title + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
    if analysis.optimiser.exitflag == 0:
        print('Saving form to:', address)
        form.to_json(address)

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
