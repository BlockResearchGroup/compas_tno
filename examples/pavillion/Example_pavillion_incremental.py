import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.algorithms import apply_sag
from compas_tno.utilities import rectangular_smoothing_constraints

# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR PAV.VAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------

t0 = thk = 0.50  # thickness on the start in meters
thk_reduction = 0.05  # in meters
thk_refined = 0.001
limit_equal = 0.005
span = 10.0  # square span for analysis
k = 1
n = 4  # Discretisation for Surfaces...
R = [5.5, 6.5, 7.5, 8.5, 9.5]
hc_list = [5.00]
# [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
sag = False
smooth = False
# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

# Basic parameters

type_structure = 'pavillionvault'
type_formdiagram = 'cross_with_diagonal'  # Try also 'fan_fd'
discretisation = 10
gradients = True

# ----------------------- Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'all',
}

form = FormDiagram.from_library(data_diagram)

if sag:
    apply_sag(form, boundary_force=sag)
if smooth:
    cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
    constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

plot_form(form, show_q=False, fix_width=10).show()

# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
# optimiser.settings['library'] = 'MMA'
# optimiser.settings['solver'] = 'MMA'
# optimiser.settings['library'] = 'IPOPT'
# optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = False
optimiser.settings['summary'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1e+3
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients

for hc in hc_list:

    print('\n**** Starting Analysis for: h={0} and thk={1} ****\n'.format(hc, thk))

    # --------------------- Shape with initial THK ---------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation*n,
        'xy_span': [[0, span], [0, k*span]],
        't': 1.0,
        'hc': hc,
        'hm': None,
        'he': None,
        'interpolation_type': 'nearest',
    }

    vault = Shape.from_library(data_shape)
    # vault.ro = 1.0

    # --------------------- Create Initial point with TNA ---------------------

    form.selfweight_from_shape(vault)
    # form = form.form_update_with_parallelisation(plot=False)
    form = form.initialise_loadpath()
    # plot_form(form, show_q=False).show()
    # from compas_plotters import MeshPlotter
    # plotter = MeshPlotter(form)
    # plotter.draw_edges()
    # plotter.draw_vertices(text={key: form.vertex_attribute(key, 'pz') for key in form.vertices()})
    # plotter.show()
    # view_thrust(form).show()
    # view_solution(form, vault).show()

    # ----------------------- Create Analysis loop on limit analysis --------------------------

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
    os.makedirs(folder, exist_ok=True)
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    if sag:
        title = title + 'sag_' + str(sag)
    if smooth:
        title = title + 'smooth_'
    forms_address = os.path.join(folder, title)

    analysis = Analysis.from_elements(vault, form, optimiser)
    results = analysis.limit_analysis_GSF(thk, thk_reduction, span, thk_refined=thk_refined, save_forms=forms_address, limit_equal=limit_equal)
    thicknesses, solutions = results

    # ----------------------- Save output data --------------------------

    csv_file = os.path.join(folder, title + '_data.csv')
    save_csv_row(thicknesses, solutions, path=csv_file, title=title)

    img_graph = os.path.join(folder, title + '_diagram.pdf')
    diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True)

    xy_limits = [[0.50, 0.15], [50, 0]]
    img_graph = os.path.join(folder, title + '_diagram_limits.pdf')
    diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits)
