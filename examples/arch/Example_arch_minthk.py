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
from copy import deepcopy
import os

# CONSTRUCT THIS!! WIP

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

H = 1.0
L = 2.0
R = H/2 + (L**2 / (8 * H))  # Note that this reduces to L/2 when springing angle = 180 deg, or L = 2*H
discretisation = 20
b = 0.5  # Out of plane dimension
type_structure = 'arch'
type_formdiagram = 'arch'

# ----------------------- 1. Create Arch shape ---------------------------

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
optimiser.data['objective'] = 'min'  # Set the objective
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1e+10  # Check if this is limiting the solution


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
