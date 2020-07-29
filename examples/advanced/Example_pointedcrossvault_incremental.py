import compas_tno
import os
import copy
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_solution
from compas_tno.plotters import diagram_of_thrust

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
thk_reduction = 0.02  # in meters
span = 10.0  # square span for analysis
R = span/2  # Only valid for rounded cross vault
k = 1/2  # Allow for rectangular Vaults

# Basic parameters

type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
discretisation = 10

# ----------------------- 1. Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
plot_form(form, show_q=False).show()

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)

# --------------------- 3. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = False
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 2000.0


# --------------------- 4. Shape with initial THK ---------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
    'hc': R*math.sqrt(1+k**2),
    'he': None,
    'hm': None,
}

vault = Shape.from_library(data_shape)
# view_shapes(vault).show()

# ----------------------- 5. Create Analysis loop on min optimisation --------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, R)
thicknesses, size_parameters, solutions_min, solutions_max = results

img_graph = os.path.join(compas_tno.get('/imgs/'), 'diagram.pdf')
diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph).show()
