import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis import Analysis
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import plot_gif_forms_and_shapes_xz
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_thrust_partial
from compas_tno.plotters import save_csv
from copy import deepcopy
import os

# ----------------------------------------------------------------------
# ---- EXAMPLE OF MIN and MAX THRUST FOR ARCH WITH INCREMENTAL THK -----
# ----------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.200  # thickness on the start in meters
thk_reduction = 0.010  # in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
xmin = []  # empty lists to keep track of  the parameters on min opt
xmax = []  # empty lists to keep track of  the parameters on max opt
forms_min = []
forms_max = []
shapes = []

while exitflag == 0:

    # Basic parameters

    H = 1.0
    L = 2.0
    R = H/2 + (L**2 / (8 * H))  # Note that this reduces to L/2 when springing angle = 180 deg, or L = 2*H
    discretisation = 20
    b = 0.5  # Out of plane dimension
    type_structure = 'arch'
    type_formdiagram = 'arch'
    t_over_R = thk/R
    print('\n----- The' , type_formdiagram, 'problem for thk/R:', t_over_R, '\n')

    # ----------------------- 1. Create Arch shape ---------------------------

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

    # ----------------------- 2. Create Form Diagram ---------------------------

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

    # --------------------- 3.1 Create Minimisation for minimum thrust ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'slsqp'
    optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
    optimiser.data['variables'] = ['ind', 'zb']
    optimiser.data['objective'] = 'min'  # Set the objective
    optimiser.data['printout'] = True
    optimiser.data['plot'] = False
    optimiser.data['find_inds'] = True
    optimiser.data['qmax'] = 1000.0  # Check if this is limiting the solution
    print('Information about optimiser:')
    print(optimiser.data, '\n')

    # --------------------------- 3.2 Run optimisation with scipy ---------------------------

    analysis = Analysis.from_elements(arch, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    # ----------------------- 4. Save and retrieve parameters ---------------------------

    file_address = compas_tno.get('test.json')
    form.to_json(file_address)
    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_min = fopt/swt
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight - min: ', fopt_over_weight_min)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        solutions_min.append(fopt_over_weight_min)
        forms_min.append(deepcopy(form))
        xmin.append(t_over_R)
        shapes.append(arch)

    # --------------------------- 5.1 Set The objective to max and run ---------------------------

    optimiser.data['objective'] = 'max'
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------------- 5.2 Extract Data from Optimisation ---------------------------

    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_max = fopt/swt
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight - max: ', fopt_over_weight_max)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        forms_max.append(deepcopy(form))
        solutions_max.append(fopt_over_weight_max)
        xmax.append(t_over_R)

    # ------------------------ 5 . Reduce the thickness ---------------------------

    if exitflag == 0:
        thk = round(thk - thk_reduction, 4)
        print('Reduce thickness to', thk, '\n')
    else:
        print('\nOptimisation did not find a solution for thk:', thk)
        print('Last solved thk:', round(thk + thk_reduction, 4))

# --------------- 6 . After exit print the list of solutions obtained ------------

print('\nSUMMARY')
print('Parameter min')
print(xmin)
print('Parameter max')
print(xmax)
print('Solutions min')
print(solutions_min)
print('Solutions max')
print(solutions_max)

# --------------- 7. Create and save the plots of the solution in data/imgs ------------


img_min = os.path.join(compas_tno.get('/imgs/'), 'min_no_voussoirs.gif')
img_max = os.path.join(compas_tno.get('/imgs/'), 'max_no_voussoirs.gif')
plot_gif_forms_and_shapes_xz(forms_min, shapes, save=img_min, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, stereotomy=discretisation).show()
plot_gif_forms_and_shapes_xz(forms_max, shapes, save=img_max, plot_reactions=True, fix_width=True, max_width=5, radius=0.02, hide_negative=True, stereotomy=discretisation).show()

img_graph = os.path.join(compas_tno.get('/imgs/'), 'diagram.pdf')
diagram_of_thrust_partial(xmin, xmax, solutions_min, solutions_max, xy_limits=[[0.20, 0.10], [60, 30]], fill=True, save=img_graph).show()

csv_path = os.path.join(compas_tno.get('/csv/'), 'diagram.csv')
save_csv(xmin, solutions_min, solutions_max, limit_state=True, path=csv_path)
