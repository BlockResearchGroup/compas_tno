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
from numpy import zeros
import math

lambd_0 = 0.05
max_lambd = 15.0
discretisation = 14
thk = 0.5
span = 10.0
spr_angle = 30.0
alpha = 1/math.cos(math.radians(spr_angle))
L = span * alpha
Ldiff = L - span
xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]

tappered = True

solve = False
visualise = True

solutions = {}

slides = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
slides = [0.8]

folder = '/Users/mricardo/compas_dev/me/hor-loads/cross/slide_diagram/minmax/'

type_diagram = 'cross_with_diagonal'

shape = Shape.create_crossvault(xy_span=xyspan_shape, thk=thk)
shape.ro = 0.1

for thk in [0.50]:

    shape_thin = Shape.create_crossvault(xy_span=xyspan_shape, thk=thk)

    for slide in slides:

        if not solve:
            continue

        # form = FormDiagram.create_parametric_form(discretisation=discretisation)
        form = FormDiagram.create_cross_form(discretisation=discretisation)
        # form = FormDiagram.create_cross_with_diagonal(discretisation=discretisation)
        # form = FormDiagram.create_fan_form()

        type_diagram = form.parameters['type']

        slide_diagram(form, delta=-slide, tappered=tappered)

        apply_selfweight_from_shape(form, shape)
        apply_horizontal_multiplier(form, lambd=lambd_0, direction='x')
        apply_bounds_on_q(form, qmax=1e-6, qmin=-20)

        # plotter = TNOPlotter(form)
        # plotter.draw_form(scale_width=False)
        # plotter.draw_supports()
        # plotter.show()

        n = form.number_of_vertices()
        load_direction = zeros((2*n, 1))
        for i, vertex in enumerate(form.vertices()):
            load_direction[i] = form.vertex_attribute(vertex, 'px')
            load_direction[i*2] = form.vertex_attribute(vertex, 'py')

        axis_sym = [[[0, 5.0], [10, 5.0]]]


        problem = Analysis.create_minthk_analysis(form, shape_thin,
        # problem = Analysis.create_minthrust_analysis(form, shape_thin,
        # problem = Analysis.create_maxthrust_analysis(form, shape_thin,
                                                    plot=False,
                                                    printout=True,
                                                    max_iter=5000,
                                                    solver='IPOPT'
                                        )
        problem.apply_envelope()
        problem.optimiser.set_features(['fixed', 'sym'])
        problem.optimiser.set_axis_symmetry(axis_sym)
        # problem.optimiser.set_additional_options(tol_inds=1e-10)
        problem.set_up_optimiser()
        # plot_svds(problem.optimiser.M)
        problem.run()

        # fopt = problem.optimiser.fopt
        thrust = form.thrust()
        swt = form.lumped_swt()

        T_over_W = abs(thrust/swt)

        form_fixed = FormDiagram.create_cross_form()
        apply_envelope_from_shape(form_fixed, shape_thin)
        shape_fixed = Shape.from_formdiagram_and_attributes(form_fixed)

        lambd_opt = abs(problem.optimiser.fopt * lambd_0)
        print('Thrust over weight:', T_over_W)

        form.attributes['T_over_W'] = T_over_W

        # view = Viewer(form, shape_fixed)
        # view.scale_edge_thickness(10.0)
        # view.draw_thrust()
        # view.draw_shape()
        # view.draw_reactions()
        # view.draw_cracks()
        # view.show()

        if problem.optimiser.exitflag == 0:
            objective = problem.optimiser.settings['objective']
            title = type_diagram + '_discr_{}_tappered_{}_slide_{}_lambd0_{}_thk_{}_{}.json'.format(discretisation, tappered, slide, lambd_0, thk, objective)
            save_json = folder + title
            form.to_json(save_json)
            print('Saved at:', save_json)
            solutions[thk] = T_over_W
        else:
            print('Problem did not find solution')

        # view = Viewer(form, shape_fixed)
        # view.scale_edge_thickness(10.0)
        # view.draw_thrust()
        # view.draw_shape()
        # view.draw_reactions()
        # view.draw_cracks()
        # view.show()


for thk in solutions:
    print(thk, solutions[thk])

# for slide in slides:

#     if not visualise:
#         continue

#     objective = problem.optimiser.settings['objective']
#     title = type_diagram + '_discr_{}_tappered_{}_slide_{}_thk_{}_{}.json'.format(discretisation, tappered, slide, thk, objective)

#     save_json = folder + title
#     form = FormDiagram.from_json(save_json)

#     form_fixed = FormDiagram.create_cross_form()
#     apply_envelope_from_shape(form_fixed, shape)
#     shape_fixed = Shape.from_formdiagram_and_attributes(form_fixed)

#     view = Viewer(form, shape_fixed)
#     view.scale_edge_thickness(10.0)
#     view.draw_thrust()
#     view.draw_shape()
#     view.draw_reactions()
#     view.draw_cracks()
#     view.show()
