import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF OPTIMISATION FOR PAVILLION VAULT MIN/MAX -------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.5
type_structure = 'pavillionvault'
type_formdiagram = 'cross_fd'
discretisation = 10

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 5.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Pavillion Vault created!')
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': discretisation,
    'fix': 'all',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
form = form.delete_boundary_edges()
print(form)
# plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 3. Create Initial point with TNA OR LP ---------------------

form = form.initialise_tna(plot=False)
# plot_form(form).show()

# optimiser = Optimiser()
# optimiser.data['library'] = 'MATLAB'
# optimiser.data['solver'] = 'SDPT3'
# optimiser.data['constraints'] = ['funicular']
# optimiser.data['variables'] = ['ind']
# optimiser.data['objective'] = 'loadpath'
# optimiser.data['printout'] = True
# optimiser.data['plot'] = False
# optimiser.data['find_inds'] = True
# optimiser.data['qmax'] = 1000.0
# analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.set_up_optimiser() # Find independent edges
# analysis.run()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 900.0
print(optimiser.data)

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

plot_form(form, show_q=False).show()

file_address = compas_tno.get('test.json')
form.to_json(file_address)

view_thrust(form).show()

# If you wish to visualise the upper and lower bound together
# view_solution(form, vault).show()
