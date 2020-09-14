import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes
import os

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF OPTIMISATION FOR PAVILLION VAULT MIN/MAX -------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.44
type_structure = 'pavillionvault'
type_formdiagram = 'cross_fd'
discretisation = 10
n=10
gradients = True

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 1.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Pavillion Vault created!')
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': discretisation*n,
    'fix': 'all',
}

form = FormDiagram.from_library(data_diagram)
# plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 3. Create Initial point with TNA OR LP ---------------------

hc = 5.0
folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
path = os.path.join(folder, title)
thk_ = 0.43
load_obj = 'max'
address_load = path + '_' + load_obj + '_thk_' + str(100*thk_) + '.json'

# form = form.initialise_tna(plot=False)
# form = form.initialise_loadpath()
form = FormDiagram.from_json(address_load)
plot_form(form, show_q=False, fix_width=False, cracks=True).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds',  'symmetry']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'max'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1e+4
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

import time
start_time = time.time()

analysis.run()

elapsed_time = time.time() - start_time
print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

file_address = path + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'
if optimiser.exitflag == 0:
    form.to_json(file_address)

plot_form(form, show_q=False, fix_width=False, cracks=True).show()
view_thrust(form).show()

# If you wish to visualise the upper and lower bound together
# view_solution(form, vault).show()
