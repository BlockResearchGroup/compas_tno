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
style_diagonals =  'straight'
t = 1.0
discretisation = [4, 12]
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
plot = True             # If selected will plot all solutions, need to be closed "by hand"

load_mult = 0.08
load_increase = 0.01
load_mult0 = load_mult
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
    't' : t           # Ammount allowed for the reactions/supports to travel down
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()

# ----------------------- 2.1. Create Form Diagram to the analysis---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': center,
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': True,
    'partial_diagonal': style_diagonals,
}

form = FormDiagram.from_library(data_diagram)
plot_form(form, show_q=False, max_width=5, simple=True).show()

# --------- 2.2. Create "Load Form Diagram" to calculate the selfweight based on the tributary area for this pattern ---------------

data_diagram = {
    'type': type_formdiagram,
    'center': center,
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}
form_load = FormDiagram.from_library(data_diagram) # This pattern is created only for the load assignment
plot_form(form_load, show_q=False).show()

# --------------------- 3. Create Initial point with LOAD PATH OPTIMISATION or TNA ---------------------

#--- If using TNA

# form = form.form_update_with_parallelisation(plot=False)
# plot_form(form).show()

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
# analysis.apply_selfweight_from_pattern(form_load, plot=True) # This adds selfweight from the "base pattern"
# analysis.set_up_optimiser()
# analysis.run()
# plot_form(form, show_q=False).show()

#--- If using load from saved form

# objective_to_load = 'min'
# load_mult_to_load = 0.04
# title = 'Dome_Px=' + str(load_mult_to_load) + '_discr_' + str(discretisation) + '_' + type_formdiagram + '_' + style_diagonals + '_' + objective_to_load
# load_json = os.path.join(folder, title + '.json')
# form = FormDiagram.from_json(load_json)
# plot_form(form, show_q=False, cracks=True, simple=True, save=os.path.join(folder, title + '.pdf')).show()

form_start = deepcopy(form)      # Keep this solution as starting point afterwards

# --------------------- 4. Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'IPOPT'
optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry-horizontal']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10e+10
print(optimiser.settings)

# --------------------------- 4.1 Set The objective to min and start loop ---------------------------

for optimiser.settings['objective'] in ['max']:

    form = form_start
    load_mult = load_mult0
    exitflag = 0

    while exitflag == 0:

        t_over_R = thk/R
        title = 'Dome_Px=' + str(load_mult) + '_discr_' + str(discretisation) + '_' + type_formdiagram + '_' + style_diagonals + '_' + optimiser.settings['objective']

        load_json = os.path.join(folder, title + '.json')
        form = FormDiagram.from_json(load_json)
        f0 = form.thrust()

        print('\n----------- Problem Title:', title, '\n')
        print('----------- F0:', f0)

        # --------------------------- 4.2. Prepare Analysis and run ---------------------------

        analysis = Analysis.from_elements(dome, form, optimiser)
        analysis.apply_selfweight_from_pattern(form_load)  # This adds selfweight from the "base pattern" does not load the intersection points
        analysis.apply_envelope()
        for key in form.vertices_where({'pz': 0.0}):  # This make sure the bounds on intrados 'lb' and extrados 'ub' do not get activated in the intersection points, i.e. these with load pz = 0.
            form.vertex_attribute(key, 'ub', radius)
            form.vertex_attribute(key, 'lb', - t)
        analysis.apply_reaction_bounds()
        if load_mult == 0.0:
            optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry']
        else:
            optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry-horizontal']
        analysis.apply_hor_multiplier(load_mult, direction_loads)  # This line applies the horizontal multiplier
        analysis.set_up_optimiser()
        analysis.run()

        # --------------------------- 4.3. Extract Data from Optimisation ---------------------------

        exitflag = optimiser.exitflag   # get info if optimisation was succeded ot not
        fopt = optimiser.fopt           # objective function optimum value
        fopt_over_weight = fopt/swt
        print('Thickness', thk)
        print('lambda:', load_mult)
        print('Thrust over weight: ', fopt_over_weight)
        print('Exitflag: ', exitflag)
        if exitflag == 0:
            # if -1 * fopt > f0:
            form.to_json(os.path.join(folder, title + '.json'))
            if optimiser.settings['objective'] == 'min':
                forms_min.append(deepcopy(form))
                solutions_min.append(fopt_over_weight)
                size_parameters_min.append(load_mult)
            else:
                forms_max.append(deepcopy(form))
                solutions_max.append(fopt_over_weight)
                size_parameters_max.append(load_mult)
            load_mult = round(load_mult + load_increase, 4)
            print('Increase Hor Multiplier to', load_mult, '\n')
            if plot:
                plot_form(form, show_q=False, cracks=True, simple=True, save=os.path.join(folder, title + '.pdf'))
        else:
            print('\nOptimisation did not find a solution for multiplier:', load_mult)
            print('Last solved multiplier:', round(load_mult - load_increase, 4))

    print('\n\nSolutions min and max so far, after running:', optimiser.settings['objective'])
    print(size_parameters_min)
    print(solutions_min)
    print(size_parameters_max)
    print(solutions_max)
