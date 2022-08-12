from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_plotters import Plotter
from numpy import zeros
import os

# Geometry parameters

radius = 5.0
thk = 0.50
solutions = {}

for np in [16]:
    solutions[np] = {}
    for nm in [20]:

        discretisation = [np, nm]
        discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

        # Parameters Optimisations

        obj = 'max_load'
        solver = 'IPOPT'
        constraints = ['funicular', 'envelope', 'reac_bounds']
        variables = ['q', 'zb', 'lambdv']
        features = ['fixed']
        starting_point = 'loadpath'
        make_video = False
        autodiff = False

        # Create shape/diagram

        dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.5)
        dome.ro = 20.0

        # form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius, diagonal=True, partial_diagonal='left')
        form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)

        # plotter = Plotter()
        # plotter.fontsize = 6
        # artist = plotter.add(form)
        # artist.draw_vertexlabels()

        # plotter.zoom_extents()
        # plotter.show()

        # form = FormDiagram.create_circular_spiral_form(discretisation=discretisation, radius=radius)

        # Maximum load magnitude

        max_load_mult = 600.0
        n = form.number_of_vertices()
        pzv = zeros((n, 1))
        pzv[0] = -1.0
        # pzv[31] = -1.0
        # pzv[101] = -1.0

        # for i in range(38):
        #     if i % 2 == 1 or i == 0:
        #         pzv[i] = -1.0

        sumpzv = float(abs(sum(pzv)))
        print('Sum of pzv:', sumpzv)

        plotter = TNOPlotter(form)
        plotter.draw_form(scale_width=False)
        plotter.draw_supports()
        plotter.show()

        # Create optimiser

        optimiser = Optimiser()
        optimiser.settings['objective'] = obj
        optimiser.settings['solver'] = solver
        optimiser.settings['constraints'] = constraints
        optimiser.settings['variables'] = variables
        optimiser.settings['features'] = features
        optimiser.settings['starting_point'] = starting_point
        optimiser.settings['derivative_test'] = False
        optimiser.settings['printout'] = False
        optimiser.settings['plot'] = False
        optimiser.settings['save_iterations'] = make_video
        optimiser.settings['autodiff'] = autodiff

        optimiser.settings['max_lambd'] = max_load_mult
        optimiser.settings['load_direction'] = pzv

        # Create analysis

        analysis = Analysis.from_elements(dome, form, optimiser)
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_reaction_bounds()

        pz0 = form.vertex_attribute(0, 'pz')

        pzt = 0
        for key in form.vertices():
            pz = form.vertex_attribute(key, 'pz')
            pzt += pz

        print('Total load of:', pzt)

        analysis.set_up_optimiser()

        # # to view starting point
        # view = Viewer(form, dome)
        # view.draw_thrust()
        # view.draw_shape()
        # view.draw_force()
        # view.show()

        analysis.run()

        fopt = analysis.optimiser.fopt
        exitflag = analysis.optimiser.exitflag

        print('Exitflag is:', exitflag)

        pc = fopt/pzt*sumpzv

        print('Percentage of load added is:', round(pc*100, 3), '%')

        if exitflag == 0:

            folder = os.path.join('/Users/mricardo/compas_dev/me/max_load/dome/apex/sensitivity', 'dome', 'radial_fd')
            os.makedirs(folder, exist_ok=True)
            title = 'dome' + '_' + 'radial_fd' + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '_pct_stw_' + str(round(pc,3)) + '.json'
            save_form = os.path.join(folder, title)
            form.to_json(save_form)

            print('Solution Saved at:', save_form)

            solutions[np][nm] = pc

    print(solutions)

# plotter = TNOPlotter(form, dome)
# plotter.show_solution()

# view = Viewer(form, dome)
# view.draw_force()
# view.show_solution()

# if make_video:

#     from compas_tno.viewers import animation_from_optimisation
#     from compas_tno.algorithms import reciprocal_from_form
#     import compas_tno

#     DATA_XFORM = compas_tno.get('Xform.json')
#     DATA_XFORCE = compas_tno.get('Xforce.json')

#     force = reciprocal_from_form(form)

#     animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, interval=150)
