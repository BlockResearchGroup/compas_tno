import os
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints


from compas.datastructures import Mesh
from compas.utilities import geometric_key

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
span = 10.0  # square span for analysis
k = 1

n = 10  # Discretisation for Surfaces...

solutions = []
hc_list = [5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]  # [5.00] # [5.00]
# [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
sag = 10
smooth = False
# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
discretisation = 10

type_topology = 'intersect'
if type_topology:
    type_formdiagram = 'topology-' + type_topology  # Try also 'fan_fd'
    discretisation = 100
    load_mesh = 'Mesh-' + type_topology
    fd_mesh = 'FormDiagram-' + type_topology

gradients = False

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
    x0 = y0 = 0
    x1 = span
    y1 = k*span
    boundary_points = [[x0, y0, 0.0], [x1, y0, 0.0], [x1, y1, 0.0], [x0, y1, 0.0]]
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
    file_ = os.path.join(folder, load_mesh + '.json')
    mesh = Mesh.from_json(file_)
    form = FormDiagram.from_mesh(mesh)
    gkey_key = form.gkey_key()
    for pt in boundary_points:
        form.vertex_attribute(gkey_key[geometric_key(pt)], 'is_fixed', True)
    form.to_json(os.path.join(folder, fd_mesh + '.json'))


if sag:
    apply_sag(form, boundary_force=sag)
if smooth:
    cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
    constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

plot_form(form, show_q=False, fix_width=10).show()

# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
# optimiser.settings['library'] = 'IPOPT'
# optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = 'target'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+10
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients

for hc in hc_list:

    print('\n**** Starting Bestfit Analysis for: h={0} and thk={1} ****\n'.format(hc, thk))

    # --------------------- Shape with initial THK ---------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation*n,
        'xy_span': [[0, span], [0, k*span]],
        't': 0.0,
        'hc': hc,
        'hm': None,
        'he': None,
    }

    vault = Shape.from_library(data_shape)

    # --------------------- Create Initial point with TNA ---------------------

    print('Calculating SWT')
    form.selfweight_from_shape(vault)
    # form = form.initialise_tna(plot=False)
    form = form.initialise_loadpath()
    # plot_form(form, show_q=False).show()

    # view_thrust(form).show()
    # view_solution(form, vault).show()

    # ----------------------- Create Analysis loop on limit analysis --------------------------

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    if load_mesh:
        title = type_structure + '_' + type_formdiagram + '_' + load_mesh + '_discr_' + str(discretisation)
    if sag:
        title = title + 'sag_' + str(sag)
    if smooth:
        title = title + 'smooth_'
    forms_address = os.path.join(folder, title)

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    solutions.append(optimiser.fopt)

    # plot_form(form, show_q=False, cracks=True).show()

    address = forms_address + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
    form.to_json(address)

print(solutions)
