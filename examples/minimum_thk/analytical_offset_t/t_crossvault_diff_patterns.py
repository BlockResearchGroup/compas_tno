
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from copy import deepcopy
import os

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

sols = {}
for discretisation in [10, 12, 14, 16, 18, 20]:

    # Basic parameters

    thk = 0.5
    span = 10.0
    k = 1.0
    n = 1
    type_structure = 'crossvault'
    type_formdiagram = 'cross_fd'
    gradients = True

    # ----------------------- 1. Create Dome shape ---------------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation*n,
        'xy_span': [[0, span], [0, k*span]],
        't': 0.0,
    }

    vault = Shape.from_library(data_shape)
    swt = vault.compute_selfweight()
    print('Selfweight computed:', swt)
    print('Vault geometry created!')

    # ----------------------- 2. Create Form Diagram ---------------------------

    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners',
    }

    form = FormDiagram.from_library(data_diagram)
    print('Form Diagram Created!')
    # print(form)
    # plot_form(form, show_q=False, fix_width=False).show()

    # --------------------- 3. Create Starting point with TNA ---------------------

    # form = form.initialise_tna(plot=False)
    form.selfweight_from_shape(vault)
    form.envelope_from_shape(vault)
    form.initialise_loadpath()
    # plot_form(form).show()

    # --------------------- 4. Create Minimisation Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'slsqp'
    optimiser.data['constraints'] = ['funicular', 'envelope']
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

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    if optimiser.exitflag == 0:
        thk_min = form.attributes['thk']
        data_shape['thk'] = thk_min
        vault = Shape.from_library(data_shape)
        form.envelope_from_shape(vault)

        folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram)
        os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
        save_form = os.path.join(folder, title)

        form.to_json(save_form + '_min_thk_' + optimiser.data['objective'] + '_' + str(thk_min) + '.json')

        sols[str(discretisation)] = thk_min

        print('Solved:', discretisation, thk_min)

    else:

        print('Not Solved:', discretisation)

print(sols)

for sol in sols:
    discr = int(sol)
    print('{0}, {1}'.format(discr, sols[sol]))
