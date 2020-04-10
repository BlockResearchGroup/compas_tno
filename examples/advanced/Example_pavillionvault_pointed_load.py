import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust_load_mult
from compas_plotters import MeshPlotter
from compas_tno.plotters import plot_forms_xz
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import plot_form_xz
from copy import deepcopy
import os

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF OPTIMISATION FOR PAVILLION VAULT MIN/MAX -------------
# ----------------------------------------------------------------------------

# Basic parameters

exitflag = 0  # means that optimisation found a solution
thk = 0.50  # thickness on the start in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
type_structure = 'pavillionvault'
type_formdiagram = 'cross_fd'
discretisation = 6
L = 10.0  # Span
R = L/2  # radius

load_mult = 2.0
load_increase = 2.0
direction_loads = 'pz'
plot_graph = True
plot_figures = True

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, L], [0.0, L]],
    't': 5.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print(type_structure, ' created!')
print('Selfweight:', swt)
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': discretisation,
    'fix': 'all',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
# form = form.delete_boundary_edges()
print(form)
# plot_form(form, show_q=False, fix_width=True).show()

# ----------------------- 1.1 Find specific node to apply load. Ex: @ Midspan or @ Quadspan ---------------------------

target_x = target_y = L/2  # Want to get the vertex having x = L/2 (midspan)
for key in form.vertices():  # Loop on all vertices
    x, y, z = form.vertex_coordinates(key)  # Get vertex coordinate
    if abs(x - target_x) < 10e-4 and abs(y - target_y) < 10e-4:  # Just to avoid any numerical error.
        target_key = key  # Assign store this key in a variable
try:
    print('Key(s) Loaded: ',target_key)
except:
    print('Did not get any key!')

plotter = MeshPlotter(form)
plotter.draw_edges()
plotter.draw_vertices(radius=0.10)
plotter.draw_vertices(keys=[target_key], radius=0.10, facecolor='FF0000')
plotter.show()  # Simple plot of the diagram

forms_min = []
forms_max = []

# --------------------- 3. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)
# plot_form(form).show()


# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
print(optimiser.data)

while exitflag == 0 and load_mult < 9.0:

    t_over_R = thk/R
    print('\n----- The' , type_formdiagram, 'problem for P=', load_mult, '\n')

    # --------------------------- 3. Prepare Analysis ---------------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.apply_pointed_load(target_key, load_mult, component=direction_loads)

    # --------------------------- 4.1 Set The objective to min and run ---------------------------

    optimiser.data['objective'] = 'min'
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------------- 4.2 Extract Data from Optimisation ---------------------------

    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_min = fopt/swt
    print('Thickness', thk)
    print('lambda:', load_mult)
    print('Thrust over weight - min: ', fopt_over_weight_min)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        forms_min.append(deepcopy(form))

    # --------------------------- 5.1 Set The objective to max and run ---------------------------

    optimiser.data['objective'] = 'max'
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------------- 5.2 Extract Data from Optimisation ---------------------------

    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_max = fopt/swt
    print('Thickness', thk)
    print('lambda:', load_mult)
    print('Thrust over weight - max: ', fopt_over_weight_max)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        forms_max.append(deepcopy(form))

    # ------------------------ 5 . Reduce the thickness or increase load ---------------------------

    if exitflag == 0:
        solutions_min.append(fopt_over_weight_min)
        solutions_max.append(fopt_over_weight_max)
        size_parameters.append(load_mult)
        load_mult = round(load_mult + load_increase, 4)
        print('Increase Hor Multiplier to', load_mult, '\n')
    else:
        print('\nOptimisation did not find a solution for multiplier:', load_mult)
        print('Last solved multiplier:', round(load_mult - load_increase, 4))


print('\n SUMMARY')
print('\nSizes calculated')
print(size_parameters)
print('Solutions Found')
print(solutions_min)
print(solutions_max)

if plot_graph:
    img_graph = os.path.join(compas_tno.get('/imgs/'),'diagram.pdf')
    diagram_of_thrust_load_mult(size_parameters, solutions_min, solutions_max, save=img_graph).show()

plot_form(form, show_q=False).show()

file_address = compas_tno.get('test_max.json')
forms_max[len(forms_max)-1].to_json(file_address)
file_address = compas_tno.get('test_min.json')
forms_min[len(forms_min)-1].to_json(file_address)

view_thrust(form).show()

# If you wish to visualise the upper and lower bound together
# view_solution(form, vault).show()
