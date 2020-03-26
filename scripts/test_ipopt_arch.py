import compas_tno
from ipopt import minimize_ipopt
# minimize_ipopt(fobj, x0, args = args, constraints = [fconstr])
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.plotters import plot_form_xz

# ----------------------------------------------------------------------
# ------ EXAMPLE OF MIN MAX THRUST FOR ARCH WITH REDUCED THK -----------
# ----------------------------------------------------------------------

# Basic parameters
# Basic parameters

thk = 0.5
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)

form = form.initialise_tna(plot=False)

# --------------------- 3.1 Create Minimisation for minimum thrust ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
print(optimiser.data)

# --------------------------- 3.2 Run optimisation with scipy ---------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

form = analysis.form
file_address = compas_tno.get('test.json')
form.to_json(file_address)
optimiser = analysis.optimiser
fopt = optimiser.fopt
view_thrust(form).show()
# plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()
