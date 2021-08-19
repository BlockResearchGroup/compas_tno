import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 2.0  # thickness on the start in meters
thk_reduction = 0.05  # in meters
solutions_min = []  # empty lists to keep track of the solutions for min thrust
solutions_max = []  # empty lists to keep track of the solutions for max thrust
size_parameters = []  # empty lists to keep track of  the parameters
thicknesses = []
span = 12.6  # square span for analysis
R = span/2  # Only valid for rounded cross vault
k = 0.6
hc = 10.0

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

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)

# --------------------- 3. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 2000.0


# --------------------- 4. Shape with initial THK ---------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
    'hc': hc,
    'hm': None,
    'he': None,
}

vault = Shape.from_library(data_shape)
# view_shapes(vault).show()

# ----------------------- 5. Create Analysis loop on limit analysis --------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, R)
thicknesses, size_parameters, solutions_min, solutions_max = results

# ----------------------- 6. Save output data --------------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', 'Amiens', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph).show()

csv_file = os.path.join(folder, title + '_data.csv')
save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)
