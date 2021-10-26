import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import Viewer
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
thk_reduction = 0.05  # in meters
thk_refined = 0.01  # thickness reduction after optimisation fails.
solutions_min = []  # empty lists to keep track of the solutions for min thrust
solutions_max = []  # empty lists to keep track of the solutions for max thrust
size_parameters = []  # empty lists to keep track of  the parameters
thicknesses = []
span = 10.0  # square span for analysis
k = 1
hc_list = [5.00]
horizontal_force_x = 10.0  # Equivalent to a
horizontal_force_y = 10.0
# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'  # Try also 'fan_fd'
discretisation = 10

# ----------------------- Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)

# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'symmetry', 'envelope', 'R']  # Note addition of constraint on rollers
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+10

for hc in hc_list:

    # --------------------- Shape with initial THK ---------------------

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

    # --------------------- Create Initial point with TNA ---------------------

    form.selfweight_from_shape(vault)
    form = form.form_update_with_parallelisation(plot=False)
    swt0 = vault.compute_selfweight()
    print('Initial Selfweight:', swt0)
    r_x = horizontal_force_x * (span * 1.0)
    r_y = horizontal_force_y * (span * 1.0)
    print('Thrust per m (x), and total:', horizontal_force_x, r_x)
    print('Thrust per m (y), and total:', horizontal_force_y, r_y)

    # ----------------------- Create Analysis loop on limit analysis --------------------------

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_partial_ratio_' + str(horizontal_force_x) + '_' + str(horizontal_force_y)
    forms_address = os.path.join(folder, title)

    analysis = Analysis.from_elements(vault, form, optimiser)
    results = analysis.limit_analysis_GSF(thk, thk_reduction, span, thk_refined=thk_refined, rollers_absolute=[r_x, r_y], save_forms=forms_address)
    thicknesses, size_parameters, solutions_min, solutions_max = results

    # ----------------------- Save output data --------------------------

    csv_file = os.path.join(folder, title + '_data.csv')
    save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)

    img_graph = os.path.join(folder, title + '_diagram.pdf')
    diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph, fill=True)
