from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import Viewer

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis
import os
from compas_tno.plotters import save_csv
from compas_tno.plotters import diagram_of_thrust

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q


radius = 5.0
radius_pattern = 3.4
# discretisation = 10
type_formdiagram = 'radial_fd'
type_structure = 'dome'
thk = 0.10
n = 6

discretisation = [4, 20]

save = False
solutions = {}

objective = ['min']
solver = 'SLSQP'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed', 'sym']
axis_sym = None  # [[[0.0, 5.0], [10.0, 5.0]], [[5.0, 0.0], [5.0, 10.0]], [[0.0, 0.0], [10.0, 10.0]]]
# qmax = 10e+6
starting_point = 'loadpath'
gradients = True

if objective == ['t']:
    variables.append(objective[0])

print(variables)

for obj in objective:  # set the objective
    solutions[obj] = {}

    for thk in [thk]:  # thickness of the problem

        # Create form diagram

        data_diagram = {
            'type': type_formdiagram,
            'center': [5.0, 5.0],
            'radius': radius_pattern,
            'discretisation': discretisation,
            'r_oculus': 0.0,
            'diagonal': False,
            'partial_diagonal': False,
        }

        form = FormDiagram.from_library(data_diagram)

        # Create shape

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': [n*discretisation[0], n*discretisation[1]],
            'center': [5.0, 5.0],
            'radius': radius,
            't': 1.0
        }

        print(data_shape)

        vault = Shape.from_library(data_shape)

        # ------------------------------------------------------------
        # -----------------------  INITIALISE   ----------------------
        # ------------------------------------------------------------

        # Apply Selfweight and Envelope

        apply_envelope_from_shape(form, vault)
        apply_selfweight_from_shape(form, vault)
        apply_bounds_on_q(form, qmax=0.0)

        form_base = form.copy()

        view = Viewer(form)
        view.view_shape()
        view.view_thrust()
        view.show()

        # ------------------------------------------------------------
        # ------------------- Proper Implementation ------------------
        # ------------------------------------------------------------

        optimiser = Optimiser()
        optimiser.settings['library'] = solver
        optimiser.settings['solver'] = solver
        optimiser.settings['constraints'] = constraints
        optimiser.settings['variables'] = variables
        optimiser.settings['features'] = features
        optimiser.settings['axis_symmetry'] = axis_sym
        optimiser.settings['objective'] = obj
        optimiser.settings['plot'] = False
        optimiser.settings['find_inds'] = False
        optimiser.settings['printout'] = True
        optimiser.settings['max_iter'] = 10000
        optimiser.settings['gradient'] = gradients
        optimiser.settings['jacobian'] = gradients
        optimiser.settings['derivative_test'] = True

        optimiser.settings['starting_point'] = starting_point

        # --------------------- 5. Set up and run analysis ---------------------

        analysis = Analysis.from_elements(vault, form, optimiser)
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_bounds_on_q()
        analysis.set_up_optimiser()
        analysis.run()

        form.overview_forces()
        if obj == 't':
            thk = form.attributes['thk']

        weight = 0
        for key in form.vertices():
            weight += form.vertex_attribute(key, 'pz')

        thrust = form.thrust()
        print('Ratio Thrust/Weight:', thrust/weight)

        address_form = '/Users/mricardo/compas_dev/compas_tno/data/form_tna_cap_min.json'
        # address_shape = '/Users/mricardo/compas_dev/compas_tno/data/shape_tna_cap.json'

        plot_form(form, show_q=False, cracks=True).show()

        view = Viewer(form)
        view.show_solution()

        form.to_json(address_form)
        # vault.to_json(address_shape)

view = Viewer(form, vault)
view.show_solution()