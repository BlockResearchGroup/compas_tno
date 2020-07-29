import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.shapes import view_shapes
from copy import deepcopy

import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.75
thk_reduction = 0.05
radius = 5.0
R = radius
type_structure = 'dome_polar'
type_formdiagram = 'radial_spaced_fd'
discretisation = [8, 20]

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*2, discretisation[1]*2],
    'center': [5.0, 5.0],
    'radius': radius,
    't' : 1.0
}

vault = Shape.from_library(data_shape)
# view_shapes(vault).show()

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

# --------------------- 3. Create Starting point with TNA ---------------------

form = form.initialise_tna(plot=False)
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1e+20
print(optimiser.data)

# ----------------------- 4. Create Analysis loop on limit analysis --------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, R)
thicknesses, size_parameters, solutions_min, solutions_max = results

# ----------------------- 5. Save output data --------------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph).show()

csv_file = os.path.join(folder, title + '_data.csv')
save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)
