import os
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import Viewer
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row
from compas_tno.algorithms import constrained_smoothing
from compas_tno.algorithms import apply_sag
from compas_tno.utilities import rectangular_smoothing_constraints

from compas_tno.problems import initialise_problem_general
from compas_tno.algorithms import z_from_form

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
type_formdiagram = 'ortho'  # Try also 'fan_fd'
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

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 1.0,
    'hc': 5.0,
    'hm': None,
    'he': None,
}

vault = Shape.from_library(data_shape)
vault.ro = 1.0

form.selfweight_from_shape(vault)

xc = 5.0
yc = 5.0
tol = 0.01

for uv in form.edges():
    p1, p2 = form.edge_coordinates(*uv)
    form.edge_attribute(uv, 'q', -1)
    if abs(p1[0] - xc) < tol and abs(p2[0] - xc) < tol:
        form.edge_attribute(uv, 'q', -5)
    elif abs(p1[1] - yc) < tol and abs(p2[1] - yc) < tol:
        form.edge_attribute(uv, 'q', -5)
    # else:
    #     pass


address = '/Users/mricardo/compas_dev/me/images/pavilion3.json'

plot_form(form, show_q=True, fix_width=10).show()

z_from_form(form)

form.to_json(address)

settings_form = {
    'color': '#D3D3D3',
    'edges.color': '#000000',
    'edges.width': 2,
    'opacity': 0.8,
    'vertices.size': 0.1,
    'vertices.on': True,
    'edges.on': True,
    'faces.on': True,
    }

# from compas_plotters import MeshPlotter
# plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
# plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
# plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
# # plotter.save(save_img)
# plotter.show()

# plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
# plotter.draw_edges(keys=[key for key in form.edges()])
# plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
# # plotter.save(save_img)
# plotter.show()

# form, force = form.reciprocal_from_form(plot=False)
# form.overview_forces()

# print('Plot of Dual')
# force.plot()

view_thrust(form, settings_form=settings_form).show()


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
