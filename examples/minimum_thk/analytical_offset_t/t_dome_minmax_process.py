
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
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [20, 16]
n = 1
ro = 1.0
gradients = True

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 0.5,
    'expanded': True
}

dome = Shape.from_library(data_shape)
dome.ro = ro
swt = dome.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')
# view_shapes(dome).show()

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
print('Form Diagram Created!')
print(form)
# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

form.envelope_from_shape(dome)
form.selfweight_from_shape(dome)
form.initialise_loadpath()
form.overview_forces()

# form = form.form_update_with_parallelisation(plot=False)
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'IPOPT'
optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb', 't']
optimiser.settings['objective'] = 't'
optimiser.settings['printout'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
print(optimiser.settings)

# --------------------- 5.1 Set up and run analysis for GSF ---------------------

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
# os.makedirs(folder, exist_ok=True)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# forms_address = os.path.join(folder, title)

# thk_max = 0.50
# thk_step = 0.05

# analysis = Analysis.from_elements(dome, form, optimiser)
# results = analysis.thk_minmax_GSF(thk_max, thk_step=thk_step, save_forms=forms_address, plot=False)
# thicknesses, solutions = results


# --------------------- 5.2 Set up and run sole analysis ---------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
os.makedirs(folder, exist_ok=True)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

forms_address = os.path.join(folder, title)

optimiser.settings['printout'] = True
# optimiser.settings['qmax'] = 5.0
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['derivative_test'] = True
optimiser.settings['objective'] = 'min'

# form = FormDiagram.from_json('/Users/mricardo/compas_dev/me/min_thk/dome/radial_fd/min_max/dome_radial_fd_discr_[20, 16]qmax=5.json')

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

print('lumped swt:', form.lumped_swt())

form.overview_forces()

from compas_plotters import MeshPlotter
plotter = MeshPlotter(form)
plotter.draw_edges()
plotter.draw_vertices(text={key: i for key, i in enumerate(form.vertices())})
plotter.show()

forms_address = os.path.join(folder, title)
form_save = forms_address + 'qmax=5' + '.json'
print(form_save)
form.to_json(form_save)
view_solution(form).show()

# ----------------------- Save output data --------------------------

# csv_file = os.path.join(folder, title + '_data.csv')
# save_csv_row(thicknesses, solutions, path=csv_file, title=title, limit_state=False)

# img_graph = os.path.join(folder, title + '_diagram.pdf')
# diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False).show()
