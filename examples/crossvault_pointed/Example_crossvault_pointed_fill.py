import compas_tno
import os
import copy
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_fill
from compas_tno.plotters import plot_form
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import save_csv

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


for fill_percentage in [1.00, 0.75]:

    exitflag = 0  # means that optimisation found a solution
    t0 = thk = 0.50  # thickness on the start in meters
    thk_reduction = 0.02  # in meters
    span = 10.0  # square span for analysis
    R = span/2  # Only valid for rounded cross vault
    k = 1.0

    # Basic parameters

    type_structure = 'pointed_crossvault'
    type_formdiagram = 'fan_fd'  # Try also 'fan_fd'
    discretisation = 10
    hc = R*math.sqrt(2)

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
    optimiser.data['qmax'] = 10e+20


    # --------------------- 4. Shape with initial THK ---------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation,
        'xy_span': [[0, span], [0, k*span]],
        't': 10.0,
        'hc': hc,
        'hm': None,
        'he': None,
    }

    vault = Shape.from_library(data_shape)
    vault.add_fill_with_height(max([point[2] for point in vault.extrados.bounding_box()]) * fill_percentage)
    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_fill_load()
    # view_shapes(vault).show()
    # view_fill(vault).show()

    # ----------------------- 5. Create Analysis loop on min optimisation --------------------------

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'SI_data', type_structure, type_formdiagram)
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_fill_' + str(fill_percentage)
    folder_forms = os.path.join(folder, 'forms', title)

    analysis = Analysis.from_elements(vault, form, optimiser)
    results = analysis.limit_analysis_GSF(thk, thk_reduction, R, fill_percentage=fill_percentage, save_forms=folder_forms, plot=False)
    thicknesses, size_parameters, solutions_min, solutions_max = results

    img_graph = os.path.join(folder, title + '_diagram.pdf')
    diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph)

    csv_file = os.path.join(folder, title + '_data.csv')
    save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)