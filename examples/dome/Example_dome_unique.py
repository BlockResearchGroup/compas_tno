import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from copy import deepcopy
import time

from compas_tno.viewers import view_shapes
import os

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.45  # 0.50
thk_ = 0.45
objective = 'max'
objective_ = 'max'
type_structure = 'dome'
type_formdiagram = 'radial_spaced_fd'
discretisation = [8, 20]
radius = 5.0
n = 10
gradients = True

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
path = os.path.join(folder, title)

# accept = False
# form = FormDiagram.from_json(os.path.join(compas_tno.get('/temp/'), 'form_' + optimiser.settings['objective'] + '.json'))
# file_address = path + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
# if accept:
#     form.to_json(file_address)
# print(x)


# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [n*discretisation[0], n*discretisation[1]],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 0.0
}


vault = Shape.from_library(data_shape)
# vault.ro = 0.1
swt = vault.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')
# view_shapes(vault).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
# optimiser.settings['library'] = 'Scipy'
# optimiser.settings['solver'] = 'slsqp'
# optimiser.settings['library'] = 'MMA'
# optimiser.settings['solver'] = 'MMA'
optimiser.settings['library'] = 'IPOPT'
optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = objective
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['summary'] = True
optimiser.settings['qmax'] = 1000.0 # swt/2
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['derivative_test'] = False

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
# print('Form Diagram Created!')
# print(form)

# --------------------- 3. Create Starting point with TNA ---------------------

if thk_:
    address_load = path + '_' + objective_ + '_thk_' + str(100*thk_) + '.json'
    form = FormDiagram.from_json(address_load)
    plot_form(form, show_q=False, fix_width=False, cracks=True).show()
    pzt = 0
    z = []
    for key in form.vertices():
        pzt += form.vertex_attribute(key, 'pz')
        z.append(form.vertex_attribute(key, 'z'))
    fac = swt/pzt
    print('Current pz: {0:.2f} | swt: {1:.2f} | Factor: {2:.2f} | z range: {3:.2f} - {4:.2f}'.format(pzt, swt, fac, min(z), max(z)))
    form.selfweight_from_shape(vault)
    form.scale_form(fac)
    pzt = 0
    z = []
    for key in form.vertices():
        pzt += form.vertex_attribute(key, 'pz')
        z.append(form.vertex_attribute(key, 'z'))
    print('Current pz: {0:.2f} | swt: {1:.2f} | Factor: {2:.2f} | z range: {3:.2f} - {4:.2f}'.format(pzt, swt, swt/pzt, min(z), max(z)))
    plot_form(form, show_q=False, fix_width=False, cracks=True).show()
else:
    form.selfweight_from_shape(vault)
    # form = form.initialise_tna(plot=False)
    form = form.initialise_loadpath()
    plot_form(form, show_q=False, fix_width=True).show()


# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

start_time = time.time()

analysis.run()

elapsed_time = time.time() - start_time
print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

form.to_json(os.path.join(compas_tno.get('/temp/'), 'form_' + optimiser.settings['objective'] + '.json'))

file_address = path + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
print('ExitFlag:', optimiser.exitflag)
if optimiser.exitflag == 0:
    form.to_json(file_address)

plot_form(form, show_q=False, fix_width=False, cracks=True).show()
# view_thrust(form).show()
