
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from compas_plotters import MeshPlotter
from copy import deepcopy
import os

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

sols = {}
for x_discr in [4, 8, 12, 16]:
    for y_discr in [12, 16, 20, 24]:
        discretisation = [x_discr, y_discr]

        thk = 0.5
        radius = 5.0
        type_structure = 'dome'
        type_formdiagram = 'radial_spaced_fd'
        # discretisation = [4, 12]
        ro = 1.0
        gradients = True
        n = 1

        # ----------------------- 1. Create Dome shape ---------------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': [discretisation[0]*n, discretisation[1]*n],
            'center': [5.0, 5.0],
            'radius': radius,
            't': 1.0
        }

        dome = Shape.from_library(data_shape)
        dome.ro = ro
        swt = dome.compute_selfweight()
        print('Selfweight computed:', swt)
        print('Vault geometry created!')

        # ----------------------- 2. Create Form Diagram ---------------------------

        data_diagram = {
            'type': type_formdiagram,
            'center': [5.0, 5.0],
            'radius': radius,
            'discretisation': discretisation,
            'r_oculus': 0.0,
            'diagonal': False,
            'partial_diagonal': False,
        }

        form = FormDiagram.from_library(data_diagram)

        # data_diagram = {
        #     'type': type_formdiagram,
        #     'center': [5.0, 5.0],
        #     'radius': radius,
        #     'discretisation': discretisation,
        #     'r_oculus': 0.0,
        #     'diagonal': False,
        #     'partial_diagonal': False,
        # }

        # form_ = FormDiagram.from_library(data_diagram)

        # plot_form(form, show_q=False, fix_width=False).show()
        # plot_form(form_, show_q=False, fix_width=False).show()

        # --------------------- 3. Create Starting point with TNA ---------------------

        # form_.selfweight_from_shape(dome)
        # form.selfweight_from_pattern(form_)
        form.envelope_from_shape(dome)
        form.selfweight_from_shape(dome)
        form.initialise_loadpath()
        # form = form.initialise_tna(plot=False)
        # plot_form(form).show()

        # --------------------- 4. Create Minimisation Optimiser ---------------------

        solvers = [['Scipy', 'SLSQP'], ['MMA', 'MMA'], ['IPOPT', 'IPOPT']]
        i = 0

        optimiser = Optimiser()
        optimiser.data['library'] = solvers[i][0]
        optimiser.data['solver'] = solvers[i][1]
        optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
        optimiser.data['variables'] = ['ind', 'zb', 't']
        optimiser.data['objective'] = 't'
        optimiser.data['printout'] = True
        optimiser.data['plot'] = False
        optimiser.data['find_inds'] = True
        optimiser.data['qmax'] = 1000.0
        optimiser.data['gradient'] = gradients
        optimiser.data['jacobian'] = gradients
        print(optimiser.data)

        # --------------------- 5. Set up and run analysis ---------------------

        analysis = Analysis.from_elements(dome, form, optimiser)
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()
        analysis.run()

        if optimiser.exitflag == 0:
            thk_min = form.attributes['thk']
            print(thk_min)
            data_shape['thk'] = thk_min
            dome = Shape.from_library(data_shape)
            form.envelope_from_shape(dome)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)

            form.to_json(save_form + '_min_thk_' + optimiser.data['objective'] + '_' + str(thk_min), '.json')

            sols[str(discretisation)] = thk_min

            print('Solved:', discretisation, thk_min)

        else:

            print('Not Solved:', discretisation)

print(sols)
