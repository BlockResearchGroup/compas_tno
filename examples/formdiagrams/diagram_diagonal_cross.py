import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.viewers.thrust import view_solution

# ----------------------------------------------------------------------------
# ---------- EXAMPLE OF MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ----------------------------------------------------------------------------

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
    'partial_bracing_modules': 4,
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()


# --------------------- 3. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)
plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'symmetry']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 10e+10
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

plot_independents(form).show()

analysis.run()

plot_form(form, show_q=False, simple=True, cracks=True).show()

file_address = compas_tno.get('test.json')
form.to_json(file_address)

view_thrust(form).show()

# If you wish to visualise the upper and lower bound together
# view_solution(form, vault).show()
