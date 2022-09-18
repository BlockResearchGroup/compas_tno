from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
# from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter

from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import slide_diagram

from compas_plotters import Plotter
from numpy import zeros
import os

# Geometry parameters

radius = 5.0
thk = 0.5
lambd0 = 0.1
max_lambd = 500.0
scale_problem = 20.0

solutions = {}

for np in [12, 16]:  # [12, 16, 20, 24]
    for nm in [12]:

        if [np, nm] in [[12, 16], [12, 20], [16, 20]]:
            continue

        discretisation = [np, nm]

        print('Analysis for:', discretisation)

        # discretisation = [12, 20]
        discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

        dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.0)
        dome.ro = scale_problem

        form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius, diagonal=True, partial_diagonal='right')
        # form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)
        type_diag = form.parameters['type']

        apply_selfweight_from_shape(form, dome)
        apply_horizontal_multiplier(form, lambd=lambd0, direction='x')

        n = form.number_of_vertices()
        load_direction = zeros((2*n, 1))
        for i, vertex in enumerate(form.vertices()):
            load_direction[i] = form.vertex_attribute(vertex, 'px')
            load_direction[i*2] = form.vertex_attribute(vertex, 'py')

        axis_sym = [[[0, 5.0], [10, 5.0]]]

        problem = Analysis.create_max_load_analysis(form, dome,
                                                    horizontal=True,
                                                    load_direction=load_direction,
                                                    max_lambd=max_lambd,
                                                    solver='IPOPT',
                                                    plot=False,
                                                    printout=True,
                                                    max_iter=5000
                                                    )
        problem.apply_envelope()
        problem.apply_reaction_bounds()
        problem.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
        problem.optimiser.set_features(['fixed'])
        problem.set_up_optimiser()

        pz0 = form.vertex_attribute(0, 'pz')

        pzt = 0
        for key in form.vertices():
            pz = form.vertex_attribute(key, 'pz')
            pzt += pz

        # print('Total load of:', pzt)


        problem.run()

        fopt = problem.optimiser.fopt
        exitflag = problem.optimiser.exitflag

        lambd = fopt * lambd0

        pc = lambd

        print('Analysis discr=', discretisation)
        print('Exitflag is:', exitflag)
        print('Optimum horizontal multiplier is:', round(pc*100, 3), '%')

        # folder = os.path.join('/Users/mricardo/compas_dev/me/max_load/dome/apex/', 'dome', 'radial_fd')
        # os.makedirs(folder, exist_ok=True)
        # title = 'dome' + '_' + 'radial_fd' + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '_pct_stw_' + str(pc) + '.json'
        # save_form = os.path.join(folder, title)
        # form.to_json(save_form)

        if exitflag == 0:
            form.attributes['lambdh'] = abs(pc)
            folder = '/Users/mricardo/compas_dev/me/hor-loads/dome/sensitivity/'
            title = 'dome_hor_load_{}_discr_{}.json'.format(type_diag, discretisation)
            save = folder + title
            form.to_json(save)
            print('Solution Saved at:', save)

            solutions[str(discretisation)] = abs(pc)

        print(solutions)

for key in solutions:
    print(key, solutions[key])
