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
from compas_tno.viewers import view_solution
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv_row

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
thk_reduction = 0.01  # in meters
solutions_min = []  # empty lists to keep track of the solutions for min thrust
solutions_max = []  # empty lists to keep track of the solutions for max thrust
size_parameters = []  # empty lists to keep track of  the parameters
thicknesses = []
span = 10.0  # square span for analysis
R = span  # Only valid for rounded cross vault
k = 1

# Basic parameters

type_structure = 'crossvault'
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
# plot_form(form, show_q=False).show()

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.form_update_with_parallelisation(plot=False)

# --------------------- 3. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'symmetry']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+10


# --------------------- 4. Shape with initial THK ---------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
    'hc': 5 * math.sqrt(2),
    'hm': None,
    'he': None
}

vault = Shape.from_library(data_shape)

# ----------------------- 5. Create Analysis loop on limit analysis --------------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
save_form = os.path.join(folder, title)

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, R, save_forms=save_form)
sizes, solutions = results

# ----------------------- 6. Save output data --------------------------


csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(sizes, solutions, path=csv_file, title=title)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(sizes, solutions, save=img_graph).show()
