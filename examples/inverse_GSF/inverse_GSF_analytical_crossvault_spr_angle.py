## Used in Computers and Structures
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from copy import deepcopy
from compas_tno.plotters import save_csv_row
from compas_tno.plotters import diagram_of_thrust
import math
import os

# ----------------------------------------------------------------------
# --------- ANALYTICAL MIN THK WITH SPRING ANGLE - CROSSVAULT ----------
# ----------------------------------------------------------------------

disc = []
discretisations = [14]
degs = [0, 20, 40]
# degs = [10]

for discretisation in discretisations:
    sols = {}

    type_formdiagram = 'fan_fd'
    type_structure = 'crossvault'
    span = span_x = span_y = 10.0
    xy_span = [[0, span_x], [0, span_y]]
    k = 1.0
    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners',
    }

    form = FormDiagram.from_library(data_diagram)
    print('Form Diagram Created!')
    j = 0

    for deg in degs:

        # Basic parameters

        # discretisation = 14
        # A = 1.00
        A = 1/math.cos(math.radians(deg))
        print('A, deg, discretisation:', A, deg, discretisation)
        xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
        thk = 0.50

        n = 4
        gradients = True

        # ----------------------- 1. Create shape ---------------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': discretisation*n,
            'xy_span': xy_span_shape,
            't': 0.0,
        }

        vault = Shape.from_library(data_shape)
        swt = vault.compute_selfweight()
        vault.ro = 1.0
        print('Selfweight computed:', swt)
        print('Vault geometry created!')

        # --------------------- 3. Create Starting point with TNA ---------------------

        # form = form.initialise_tna(plot=False)
        form.selfweight_from_shape(vault)
        # if j == 0:
        form.initialise_loadpath()
        # plot_form(form).show()

        # --------------------- 4. Create Minimisation Optimiser ---------------------

        optimiser = Optimiser()
        optimiser.data['library'] = 'Scipy'
        optimiser.data['solver'] = 'slsqp'
        optimiser.data['constraints'] = ['funicular', 'envelope']
        optimiser.data['variables'] = ['ind', 'zb', 't']
        optimiser.data['objective'] = 't'
        optimiser.data['min_thk'] = 0.0
        optimiser.data['max_thk'] = thk*1.0
        optimiser.data['printout'] = False
        optimiser.data['plot'] = False
        optimiser.data['find_inds'] = True
        optimiser.data['qmax'] = 1000.0
        optimiser.data['gradient'] = gradients
        optimiser.data['jacobian'] = gradients

        # --------------------- 5. Set up and run analysis ---------------------

        folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'deg='+str(deg))
        os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_deg=' + str(deg)
        forms_address = os.path.join(folder, title)

        thk_max = thk
        thk_step = 0.05

        analysis = Analysis.from_elements(vault, form, optimiser)
        results = analysis.thk_minmax_GSF(thk_max, thk_step=thk_step, save_forms=forms_address, plot=False)
        thicknesses, solutions = results

        # ----------------------- Save output data --------------------------

        csv_file = os.path.join(folder, title + '_data.csv')
        save_csv_row(thicknesses, solutions, path=csv_file, title=title, limit_state=False)

        img_graph = os.path.join(folder, title + '_diagram.pdf')
        diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False).show()

        # analysis = Analysis.from_elements(vault, form, optimiser)
        # analysis.apply_selfweight()
        # analysis.apply_envelope()
        # analysis.apply_reaction_bounds()
        # analysis.set_up_optimiser()
        # analysis.run()

        j += 1

        if optimiser.exitflag == 0:
            thk_min = form.attributes['thk']
            data_shape['thk'] = thk_min
            vault = Shape.from_library(data_shape)
            form.envelope_from_shape(vault)

            # folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram)
            # os.makedirs(folder, exist_ok=True)
            # title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_deg=' + str(deg)
            # save_form = os.path.join(folder, title)

            # form.to_json(save_form + '_min_thk_' + optimiser.data['objective'] + '_' + str(thk_min) + '.json')

            # sols[str(discretisation)] = thk_min
            sols[str(deg)] = thk_min
            # sols[str(A)] = thk_min


            # print('Solved:', discretisation, thk_min)
            print('Solved:', deg, thk_min)

        else:

            # print('Not Solved:', discretisation)
            print('Not Solved:', deg)

        # plot_form(form, show_q=False, simple=True, cracks=True).show()
        # view_solution(form, vault).show()

    print(discretisation)
    print(sols)
    disc.append(sols)

print(disc)

i = 0
for sols in disc:
    print('discretisation:', discretisations[i])
    for sol in sols:
        discr = float(sol)
        print('{0}, {1}'.format(discr, sols[sol]))
    i += 1
