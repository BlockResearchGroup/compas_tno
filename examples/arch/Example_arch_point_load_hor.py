import compas_tno
import os
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis.analysis import Analysis
from compas_tno.plotters import diagram_of_thrust_load_mult
from compas_plotters import MeshPlotter
from compas_tno.plotters import plot_forms_xz
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import plot_form_xz
from copy import deepcopy

# ----------------------------------------------------------------------
# ---- EXAMPLE OF MIN and MAX THRUST FOR ARCH WITH HIGHER POINT LOAD -----
# ----------------------------------------------------------------------

# Setup Parameters

exitflag = 0  # means that optimisation found a solution
thk = 0.200  # thickness on the start in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
load_mult = 0.0
load_increase = 0.5
direction_loads = 'px'
plot_graph = True
plot_figures = True


# Basic parameters

H = 1.0
L = 2.0
R = H/2 + (L**2 / (8 * H))  # Note that this reduces to L/2 when springing angle = 180 deg, or L = 2*H
discretisation = 25  # Uneven number of nodes, so it exists a node in the midspan.
b = 0.5  # Out of plane dimension
type_structure = 'arch'
type_formdiagram = 'arch'

# ----------------------- 1. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'H': H,
    'L': L,
    'total_nodes': discretisation,
    'x0': 0.0
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created, see overview:')
print(form)

# ----------------------- 1.1 Find specific node to apply load. Ex: @ Midspan or @ Quadspan ---------------------------

target_x = L/2  # Want to get the vertex having x = L/2 (midspan)
for key in form.vertices():  # Loop on all vertices
    x, y, z = form.vertex_coordinates(key)  # Get vertex coordinate
    if abs(x - target_x) < 10e-4:  # Just to avoid any numerical error.
        target_key = key  # Assign store this key in a variable
try:
    print('Key(s) Loaded: ',target_key)
except:
    print('Did not get any key!')

plotter = MeshPlotter(form)
plotter.draw_edges()
plotter.draw_vertices(radius=0.01)
plotter.draw_vertices(keys=[target_key], radius=0.01, facecolor='FF0000')
plotter.show()  # Simple plot of the diagram

forms_min = []
forms_max = []

# --------------------- 2. Create Optimiser still with no objective  ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 5000.0  # Check if this is limiting the solution

while exitflag == 0:

    t_over_R = thk/R
    print('\n----- The' , type_formdiagram, 'problem for thk/R:', t_over_R, 'and P=', load_mult, '\n')

    # ----------------------- 1. Create Arch shape for thk ---------------------------

    data_shape = {
        'type': type_structure,
        'H': H,
        'L': L,
        'thk': thk,
        'discretisation': discretisation,
        'b': b,
        't': 10.0,
        'x0': 0.0
    }

    arch = Shape.from_library(data_shape)
    area = arch.middle.area()
    swt = arch.compute_selfweight()
    print('Arch created!')
    print('Self-weight is:', swt)
    print('Area is:', area)
    # view_shapes(arch).show()

    # --------------------------- 3. Prepare Analysis ---------------------------

    analysis = Analysis.from_elements(arch, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.apply_pointed_load(target_key, -1*load_mult*form.vertex_attribute(target_key, 'pz'), proportional=False, component=direction_loads)

    # --------------------------- 4.1 Set The objective to min and run ---------------------------

    optimiser.settings['objective'] = 'min'
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
    # plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()

    # --------------------------- 5.1 Set The objective to max and run ---------------------------

    optimiser.settings['objective'] = 'max'
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
    # plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()

    # ------------------------ 5 . Reduce the thickness ---------------------------

    if exitflag == 0:
        solutions_min.append(fopt_over_weight_min)
        solutions_max.append(fopt_over_weight_max)
        size_parameters.append(load_mult)
        load_mult = round(load_mult + load_increase, 4)
        print('Increase Hor Multiplier to', load_mult, '\n')
    else:
        print('\nOptimisation did not find a solution for multiplier:', load_mult)
        print('Last solved multiplier:', round(load_mult - load_increase, 4))

# --------------- 6 . After exit print the list of solutions obtained ------------

img_file_min = os.path.join(compas_tno.get('/imgs/'),'test_min.gif')
img_file_max = os.path.join(compas_tno.get('/imgs/'),'test_max.gif')
# plot_gif_forms_xz(forms_min, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, save=img_file_min).show()
# plot_gif_forms_xz(forms_max, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, save=img_file_max).show()
# plot_forms_xz(forms_min, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()
# plot_forms_xz(forms_max, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()
print('\n SUMMARY')
print('\nSizes calculated')
print(size_parameters)
print('Solutions Found')
print(solutions_min)
print(solutions_max)

if plot_figures:
    img_file_min = os.path.join(compas_tno.get('/imgs/'),'test_min.gif')
    img_file_max = os.path.join(compas_tno.get('/imgs/'),'test_max.gif')
    plot_gif_forms_xz(forms_min, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, save=img_file_min).show()
    plot_gif_forms_xz(forms_max, arch, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, save=img_file_max).show()
if plot_graph:
    img_graph = os.path.join(compas_tno.get('/imgs/'),'diagram.pdf')
    diagram_of_thrust_load_mult(size_parameters, solutions_min, solutions_max, save=img_graph).show()
