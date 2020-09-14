import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.algorithms import apply_sag
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.plotters import open_csv_row
from compas.datastructures import Mesh
from compas.utilities import geometric_key
from compas_tno.algorithms import constrained_smoothing
from compas_tno.utilities import rectangular_smoothing_constraints
from compas_triangle.delaunay import conforming_delaunay_triangulation


# ------------------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ------------------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.60  # thickness on the start in meters
thk_reduction = 0.05  # in meters
thk_refined = 0.01
solutions_min = []  # empty lists to keep track of the solutions for min thrust
solutions_max = []  # empty lists to keep track of the solutions for max thrust
size_parameters = []  # empty lists to keep track of  the parameters
thicknesses = []
span = 10.0  # square span for analysis
k = 1
hc_list = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]  # [5.00] # [5.00]
sag = False
smooth = False
# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]
plot = False
gradients = True

# Basic parameters

type_structure = 'pointed_crossvault'
type_topology = 'domical'
type_formdiagram = 'topology-' + type_topology  # Try also 'fan_fd'
discretisation = 100
load_mesh = 'Mesh-' + type_topology
starting_point = 'loadpath'

# ----------------------- Create Form Diagram for analysis ---------------------------

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

# if load_mesh == 'triangle':
#     form = FormDiagram.from_triangle(boundary_points)
#     gkey_key = form.gkey_key()

if sag:
    apply_sag(form, boundary_force=sag)

if smooth:
    cons = rectangular_smoothing_constraints(form, xy_span=[[0, span], [0, k*span]])
    constrained_smoothing(form, damping=0.5, kmax=100, constraints=cons, algorithm='centroid')

form.to_json(os.path.join(folder, 'FormDiagram-' + type_topology + '.json'))

plot_form(form, show_q=False, fix_width=10).show()


# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 10e+10
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

for hc in hc_list:

    # --------------------- Shape with initial THK ---------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation,
        'xy_span': [[0, span], [0, k*span]],
        't': 0.0,
        'hc': hc,
        'hm': None,
        'he': None,
    }

    vault = Shape.from_library(data_shape)

    # ----------------------- Create title for analysis --------------------------

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
    title = type_structure + '_' + type_formdiagram + '_' + load_mesh + '_discr_' + str(discretisation)
    if sag:
        title = title + 'sag_' + str(sag)
    if smooth:
        title = title + 'smooth_'
    forms_address = os.path.join(folder, title)

    # --------------------- Create Starting Point ---------------------

    form.selfweight_from_shape(vault)

    if starting_point == 'tna':
        form = form.initialise_tna(plot=False)
    elif starting_point == 'loadpath':
        form = form.initialise_loadpath()
    else:
        csv_file = os.path.join(folder, title + '_data.csv')
        sizes0, solutions0 = open_csv_row(csv_file, printout=False)
        thk = sizes0[0][-1]
        print('------ Starting from thickness t=', thk)
        thk_reduction = 0.01
        thk_refined = 0.001
        address = forms_address + '_' + 'min' + '_thk_' + str(100*thk) + '.json'
        form = FormDiagram.from_json(address)

    plot_form(form, show_q=False).show()
    # view_thrust(form).show()

    # ----------------------- Create Analysis loop on limit analysis --------------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    results = analysis.limit_analysis_GSF(thk, thk_reduction, span, thk_refined=thk_refined, save_forms=forms_address, plot=plot)
    sizes, solutions = results

    # ----------------------- Save output data --------------------------

    csv_file = os.path.join(folder, title + '_data.csv')
    # save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)
    # sizes = [sizes0[0] + sizes[0][1:], sizes0[1] + sizes[1][1:]]
    # solutions = [solutions0[0] + solutions[0][1:], solutions0[1] + solutions[1][1:]]
    save_csv_row(sizes, solutions, path=csv_file, title=title)

    xy_limits = [[0.60, 0.20], [120, 30]]
    img_graph = os.path.join(folder, title + '_diagram.pdf')
    diagram_of_thrust(sizes, solutions, save=img_graph, fill=True, xy_limits=xy_limits).show()
