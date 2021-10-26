import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import Viewer
from compas_tno.viewers import view_shapes
from compas_tno.plotters import diagram_of_thrust_load_mult
from compas_plotters import MeshPlotter
from compas_tno.plotters import plot_forms_xz
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import plot_form_xz
from copy import deepcopy
import os

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF OPTIMISATION FOR DOME WITH HOR. MULTIPLIER -------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.50
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 20]
center = [5.0, 5.0]
R = radius

exitflag = 0        # means that optimisation found a solution
thk = 0.50             # thickness on the start in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
forms_min = []
forms_max = []
size_parameters_min = []  # empty lists to keep track of the parameters
size_parameters_max = []  # empty lists to keep track of the parameters
size_parameters = []

load_mult = 0.20
load_increase = 0.01
direction_loads = 'px'
plot_graph = True
plot_figures = True
folder = compas_tno.get('/dome/') # Folder to Save the structure

print(folder)

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*2, discretisation[1]*2],
    'center': center,
    'radius': radius,
    't' : 1.0           # Ammount allowed for the reactions/supports to travel down
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()

# ----------------------- 2. Create Form Diagram and "Load Form Diagram" ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': center,
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': True,
    'partial_diagonal': 'left',
}

form = FormDiagram.from_library(data_diagram)
plot_form(form, show_q=False).show()

# --------------------- 3. Create Initial point with LOAD PATH OPTIMISATION or TNA ---------------------

#--- If using TNA

form = form.form_update_with_parallelisation(plot=False)
plot_form(form).show()

#--- If using LOADPATH

# optimiser = Optimiser()
# optimiser.settings['library'] = 'MATLAB'
# optimiser.settings['solver'] = 'SDPT3'
# optimiser.settings['constraints'] = ['funicular']
# optimiser.settings['variables'] = ['ind']
# optimiser.settings['objective'] = 'loadpath'
# optimiser.settings['printout'] = True
# optimiser.settings['plot'] = False
# optimiser.settings['find_inds'] = True
# optimiser.settings['qmax'] = 10e+10
# analysis = Analysis.from_elements(dome, form, optimiser)
# analysis.apply_selfweight()
# analysis.set_up_optimiser()
# analysis.run()

SWT = dome.compute_selfweight()
PZSUM = 0
for key in form.vertices():
    PZSUM += form.vertex_attribute(key, 'pz')
print(SWT, PZSUM)

plot_form(form, show_q=False).show()


#--- If using load from saved form

# title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + 'min'
# load_json = os.path.join(folder, title + '.json')
# form = FormDiagram.from_json(load_json)

# plot_form(form, show_q=False).show()

# form_start = deepcopy(form)      # Keep this solution stored to perhaps use it as starting point afterwards

# --------------------- 4. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds' 'symmetry-horizontal']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+4
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
print(optimiser.settings)

# --------------------------- 4.1 Set The objective to min and start loop ---------------------------

optimiser.settings['objective'] = 'min'

while exitflag == 0:

    t_over_R = thk/R
    title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective']

    print('\n----------- Problem Title:', title, '\n')

    # --------------------------- 4.2. Prepare Analysis and run ---------------------------

    analysis = Analysis.from_elements(dome, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.apply_hor_multiplier(load_mult, direction_loads)  # This line applies the horizontal multiplier
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------------- 4.3. Extract Data from Optimisation ---------------------------

    exitflag = optimiser.exitflag   # get info if optimisation was succeded ot not
    fopt = optimiser.fopt           # objective function optimum value
    fopt_over_weight_min = fopt/swt
    print('Thickness', thk)
    print('lambda:', load_mult)
    print('Thrust over weight - min: ', fopt_over_weight_min)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        form.to_json(os.path.join(folder, title + '.json'))
        forms_min.append(deepcopy(form))
        solutions_min.append(fopt_over_weight_min)
        size_parameters_min.append(load_mult)
        load_mult = round(load_mult + load_increase, 4)
        print('Increase Hor Multiplier to', load_mult, '\n')
    else:
        print('\nOptimisation did not find a solution for multiplier:', load_mult)
        print('Last solved multiplier:', round(load_mult - load_increase, 4))


#--------------------------- 5.0 If want to start from the same starting point ---------------------------

# form = form_start

# --------------------------- 5.2 Set The objective to max and run ---------------------------

optimiser.settings['objective'] = 'max'
exitflag = 0
load_mult = 0.0           # Start from zero Horizontal Load applied, or maybe from a certain amount...

while exitflag == 0:

    t_over_R = thk/R
    title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective']

    print('\n----------- Problem Title:', title, '\n')

    # --------------------------- 4.2. Prepare Analysis ---------------------------

    analysis = Analysis.from_elements(dome, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.apply_hor_multiplier(load_mult, direction_loads)
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------------- 4.2 Extract Data from Optimisation ---------------------------

    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_max = fopt/swt
    print('Thickness', thk)
    print('lambda:', load_mult)
    print('Thrust over weight - max: ', fopt_over_weight_max)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        form.to_json(os.path.join(folder, title + '.json'))
        forms_max.append(deepcopy(form))
        solutions_max.append(fopt_over_weight_max)
        size_parameters_max.append(load_mult)
        load_mult = round(load_mult + load_increase, 4)
        print('Increase Hor Multiplier to', load_mult, '\n')
    else:
        print('\nOptimisation did not find a solution for multiplier:', load_mult)
        print('Last solved multiplier:', round(load_mult - load_increase, 4))


print('\n SUMMARY')
print('\nSizes calculated')
print(size_parameters_min)
print(size_parameters_max)
print('Solutions Found')
print(solutions_min)
print(solutions_max)
