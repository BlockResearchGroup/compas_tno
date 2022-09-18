from compas_tno import analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.utilities import slide_diagram
from compas_tno.analysis import Analysis
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.problems.problems import plot_svds
from compas_tno.viewers import Viewer
from compas_tno.algorithms import apply_sag
from numpy import zeros
import math

lambd_0 = 1.0
max_lambd = 10.0
discretisation = 14
span = 10.0
spr_angle = 30.0
alpha = 1/math.cos(math.radians(spr_angle))
L = span * alpha
Ldiff = L - span
xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]

tappered = True

solve = True
visualise = True

solutions = {}

slides = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
slides = [0.5]

shape = Shape.create_crossvault(xy_span=xyspan_shape)
shape.ro = 0.1

folder = '/Users/mricardo/compas_dev/me/hor-loads/cross/slide_diagram/'

type_diagram = 'cross_with_diagonal'

for slide in slides:

    if not solve:
        continue

    # form = FormDiagram.create_parametric_form(discretisation=discretisation)
    form = FormDiagram.create_cross_form(discretisation=discretisation)
    # form = FormDiagram.create_cross_with_diagonal(discretisation=discretisation)
    # form = FormDiagram.create_fan_form()

    # type_diagram = form.parameters['type']
    # type_diagram = 'mixed'
    # type_diagram = 'sag'

    # path = '/Users/mricardo/compas_dev/me/hor-loads/cross/mixed_form3.json'
    # form = FormDiagram.from_json(path)

    # move_pattern_to_origin(form)

    # title = 'cross_discr_{}_full_slide_{}_lambdh.json'.format(discretisation, slide)
    # title = type_diagram + '_discr_{}_tappered_{}_slide_{}_lambdh.json'.format(discretisation, tappered, slide)

    slide_diagram(form, delta=-slide, tappered=tappered)

    apply_selfweight_from_shape(form, shape)
    apply_horizontal_multiplier(form, lambd=lambd_0, direction='x')
    apply_bounds_on_q(form, qmax=1e-6, qmin=-20)

    # apply_sag(form, boundary_force=10)

    # apply_horizontal_multiplier(form, lambd=0.1, direction='x')

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_supports()
    plotter.show()

    n = form.number_of_vertices()
    load_direction = zeros((2*n, 1))
    for i, vertex in enumerate(form.vertices()):
        load_direction[i] = form.vertex_attribute(vertex, 'px')
        load_direction[i*2] = form.vertex_attribute(vertex, 'py')

    axis_sym = [[[0, 5.0], [10, 5.0]]]

    problem = Analysis.create_max_load_analysis(form, shape,
                                                horizontal=True,
                                                load_direction=load_direction,
                                                max_lambd=max_lambd,
                                                solver='IPOPT',
                                                plot=True,
                                                printout=True,
                                                max_iter=5000
                                                )
    problem.apply_envelope()
    problem.optimiser.set_features(['fixed', 'sym'])
    problem.optimiser.set_axis_symmetry(axis_sym)
    # problem.optimiser.set_additional_options(tol_inds=1e-6)
    problem.optimiser.set_additional_options(solver_convex='MATLAB')
    problem.set_up_optimiser()
    plot_svds(problem.optimiser.M)
    problem.run()

    form_fixed = FormDiagram.create_cross_form()
    apply_envelope_from_shape(form_fixed, shape)
    shape_fixed = Shape.from_formdiagram_and_attributes(form_fixed)

    lambd_opt = abs(problem.optimiser.fopt * lambd_0)
    print('Horizontal_mulltiplier:', lambd_opt)

    form.attributes['lambdh'] = lambd_opt

    # view = Viewer(form, shape_fixed)
    # view.scale_edge_thickness(10.0)
    # view.draw_thrust()
    # view.draw_shape()
    # view.draw_reactions()
    # view.draw_cracks()
    # view.show()

    # if problem.optimiser.exitflag == 0:
    #     save_json = folder + title
    #     form.to_json(save_json)

    #     print('Saved at:', save_json)

    #     solutions[slide] = lambd_opt
    # else:
    #     print('Problem did not find solution')

# # for slide in solutions:
# #     print(slide, solutions[slide])

# for slide in slides:

#     if not visualise:
#         continue

#     type_diagram = 'mixed'
#     title = type_diagram + '_discr_{}_tappered_{}_slide_{}_lambdh.json'.format(discretisation, tappered, slide)

#     save_json = folder + title
#     form = FormDiagram.from_json(save_json)

#     form_fixed = FormDiagram.create_cross_form()
#     apply_envelope_from_shape(form_fixed, shape)
#     shape_fixed = Shape.from_formdiagram_and_attributes(form_fixed)

    view = Viewer(form, shape_fixed)
    view.scale_edge_thickness(10.0)
    view.draw_thrust()
    view.draw_shape()
    view.draw_reactions()
    view.draw_cracks()
    view.show()
