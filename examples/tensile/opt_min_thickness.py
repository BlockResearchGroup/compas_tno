from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.viewers import view_solution

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 10 * discretisation

save = False
solutions = {}

objective = ['t']  # try 'max'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 't']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'

for obj in objective:  # set the objective

    # Create form diagram

    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners'
    }

    form = FormDiagram.from_library(data_diagram)

    # Create shape

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation_shape,
        'xy_span': [[0, span], [0, k*span]],
        'center': [5.0, 5.0],
        'radius': span/2,
        't': 0.0,
    }

    vault = Shape.from_library(data_shape)
    vault.ro = 20.0

    # ------------------------------------------------------------
    # -----------------------  INITIALISE   ----------------------
    # ------------------------------------------------------------

    # Apply Selfweight and Envelope

    apply_envelope_from_shape(form, vault)
    apply_selfweight_from_shape(form, vault)
    apply_bounds_on_q(form, qmax=1e-6)

    # ------------------------------------------------------------
    # ------------------- Proper Implementation ------------------
    # ------------------------------------------------------------

    optimiser = Optimiser()
    optimiser.settings['library'] = solver
    optimiser.settings['solver'] = solver
    optimiser.settings['constraints'] = constraints
    optimiser.settings['variables'] = variables
    optimiser.settings['features'] = features
    optimiser.settings['objective'] = obj
    optimiser.settings['plot'] = False
    optimiser.settings['find_inds'] = False
    optimiser.settings['max_iter'] = 500
    optimiser.settings['gradient'] = True
    optimiser.settings['jacobian'] = True
    optimiser.settings['printout'] = True
    optimiser.settings['starting_point'] = starting_point
    optimiser.settings['sym_loads'] = False

    # --------------------- 5. Set up and run analysis ---------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.set_up_optimiser()
    analysis.run()

    plot_form(form).show()
    view_solution(form).show()
