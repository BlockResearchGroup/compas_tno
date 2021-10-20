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
from compas_tno.viewers import Viewer
from compas_tno.viewers import  view_shapes
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row

from compas_tno.algorithms import constrained_smoothing
from compas_tno.algorithms import apply_sag
from compas_tno.utilities import rectangular_smoothing_constraints

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.80  # thickness on the start in meters
span = 10.0  # square span for analysis
k = 1

n = 5  # Discretisation for Surfaces...

thk_reduction = 0.05  # in meters
thk_refined = 0.0005
limit_equal = 0.001

solutions = []
sag = False
smooth = True

# Choice of Form Diagram

type_structure = 'domicalvault'
type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
discretisation = 10
gradients = False

type_topology = None  # None
if type_topology:
    type_formdiagram = 'topology-' + type_topology
    discretisation = 50
    fd_mesh = 'FormDiagram-' + type_topology

# ----------------------- Create Form Diagram for analysis ---------------------------

if not type_topology:
    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners',
    }
    form = FormDiagram.from_library(data_diagram)
else:
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
    form = FormDiagram.from_json((os.path.join(folder, fd_mesh + '.json')))

if sag:
    apply_sag(form, boundary_force=sag)
if smooth:
    cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
    constrained_smoothing(form, damping=0.5, kmax=500, constraints=cons, algorithm='centroid')

plot_form(form, show_q=False, fix_width=10).show()

# --------------------- Shape with initial THK and SWR ---------------------

data = {
    'type': type_structure,
    'xy_span': [[0, span], [0, k*span]],
    'thk': thk,
    'discretisation': [n*discretisation]*2,
    # 'center': [5.0, 5.0],
    # 'radius': 5.0,
    't': 0.0,
}

vault = Shape.from_library(data)
form.selfweight_from_shape(vault)
print('SWT Calculated')

# --------------------- Create Initial point with TNA ---------------------

# form = form.form_update_with_parallelisation(plot=False)
form = form.initialise_loadpath()
# view_thrust(form).show()

# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+10
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients


# If bestfit desirable:
# optimiser.settings['objective'] = 'target'
# analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()
# analysis.run()
# f = form.distance_target()
# print('Best-fit Distance -', type_formdiagram, ':', f)

# ----------------------- 5. Create Analysis loop on limit analysis --------------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
forms_address = os.path.join(folder, title)

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.limit_analysis_GSF(thk, thk_reduction, span, thk_refined=thk_refined, save_forms=forms_address, limit_equal=limit_equal)
thicknesses, solutions = results
print(results)

# ----------------------- 6. Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(thicknesses, solutions, path=csv_file, title=title)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True)

xy_limits = [[0.60, 0.20], [120, 30]]
img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits)
