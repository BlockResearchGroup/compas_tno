
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from copy import deepcopy
import os
from compas_tno.plotters import save_csv_row
from compas_tno.plotters import diagram_of_thrust


# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
span = 10.0
k = 1.0
n = 1
# ro =  1.0
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
gradients = True
sag = None
smooth = None

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

vault = Shape.from_library(data_shape)
# vault.ro = ro
swt = vault.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb', 't']
optimiser.data['objective'] = 't'
optimiser.data['printout'] = False
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

hc = 5.0
folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
os.makedirs(folder, exist_ok=True)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
if sag:
    title = title + 'sag_' + str(sag)
if smooth:
    title = title + 'smooth_'
forms_address = os.path.join(folder, title)

thk_max = 0.50
thk_step = 0.05

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.thk_minmax_GSF(thk_max, thk_step=thk_step, save_forms=forms_address, plot=False)
thicknesses, solutions = results

# ----------------------- Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(thicknesses, solutions, path=csv_file, title=title, limit_state=False)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False).show()
