import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.shapes import view_shapes
from copy import deepcopy

from compas_tno.algorithms import apply_sag
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints

import os
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import open_csv
from compas_tno.plotters import save_csv_row

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.75
thk_reduction = 0.05  # in meters
thk_refined = 0.001
limit_equal = 0.005
radius = 5.0
span = radius * 2
R = radius
n = 10  # Discretisation for Surfaces...
sag = False
smooth = False

type_structure = 'dome'
type_formdiagram = 'radial_spaced_fd'
discretisation = [8, 20]
gradients = True

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 1.0,
    'interpolation_type': 'linear',
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

if sag:
    apply_sag(form, boundary_force=sag)
if smooth:
    cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
    constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

plot_form(form, show_q=False, fix_width=10).show()

# --------------------- 3. Create Starting point with TNA ---------------------

form.selfweight_from_shape(vault)
# form = form.form_update_with_parallelisation(plot=False)
form = form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = False
optimiser.settings['summary'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1e+8
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['summary'] = True

# ----------------------- Create Analysis loop on limit analysis --------------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
os.makedirs(folder, exist_ok=True)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
if sag:
    title = title + 'sag_' + str(sag)
if smooth:
    title = title + 'smooth_'
forms_address = os.path.join(folder, title)

# ----------------------- 4. Create Analysis loop on limit analysis --------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, span, thk_refined=thk_refined, save_forms=forms_address, limit_equal=limit_equal)
thicknesses, solutions = results

# ----------------------- 5. Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(thicknesses, solutions, path=csv_file, title=title)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True).show()

xy_limits = [[0.50, 0.15], [100, 0]]
img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits).show()
