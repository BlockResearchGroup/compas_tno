import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis import Analysis
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_intrados
from copy import deepcopy
import os
from compas.datastructures import Mesh

# ----------------------------------------------------------------------
# ---- EXAMPLE OF MIN and MAX THRUST FOR ARCH WITH INCREMENTAL THK -----
# ----------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.20  # thickness on the start in meters
thk_reduction = 0.010  # in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
forms_min = []
forms_max = []

test_shape = 1

H = 1.0
L = 2.0
R = H/2 + (L**2 / (8 * H))  # Note that this reduces to L/2 when springing angle = 180 deg, or L = 2*H
discretisation = 20
b = 0.5  # Out of plane dimension
type_structure = 'arch'
type_formdiagram = 'arch'

# ----------------------- 1. Create Arch shape ---------------------------

folder_input = compas_tno.get('/input/')
print(folder_input)

intrados = Mesh.from_json(os.path.join(folder_input, 'intrados{0}.json'.format(test_shape)))
extrados = Mesh.from_json(os.path.join(folder_input, 'extrados{0}.json'.format(test_shape)))
middle = Mesh.from_json(os.path.join(folder_input, 'target.json'.format(test_shape)))

data_shape = {
    'type': type_structure,
    'H': H,
    'L': L,
    'thk': thk,
    'discretisation': discretisation,
    'b': b,
    't': 5.0,
    'x0': 0.0
}

vault = Shape.from_library(data_shape)
vault = Shape.from_meshes(intrados, extrados, middle, data=data_shape)
swt = vault.compute_selfweight()
print('selfweight:', swt)
# view_shapes(vault).show()


# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'H': H,
    'L': L,
    'total_nodes': discretisation,
    'x0': 0.0
}

form = FormDiagram.from_library(data_diagram)

# --------------------- 3 Create Minimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'max'  # Set the objective
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1e+10  # Check if this is limiting the solution


analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()
