import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis import Analysis
from compas_tno.plotters import plot_gif_forms_xz
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_intrados
from copy import deepcopy
import os
from compas.datastructures import Mesh

# ----------------------------------------------------------------------
# ---- EXAMPLE OF MIN and MAX THRUST FOR ARCH WITH INCREMENTAL THK -----
# ----------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.20  # thickness on the start in meters
thk_reduction = 0.010  # in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
forms_min = []
forms_max = []
shapes = []

test_shape = 2

# Get the meshes on the bounds

folder_input = compas_tno.get('/input/')
print(folder_input)

intrados = Mesh.from_json(os.path.join(folder_input, 'intrados{0}.json'.format(test_shape)))
extrados = Mesh.from_json(os.path.join(folder_input, 'extrados{0}.json'.format(test_shape)))
middle = Mesh.from_json(os.path.join(folder_input, 'target.json'.format(test_shape)))

while exitflag == 0:

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
        't': 5.0,
        'x0': 0.0
    }

    vault = Shape.from_library(data_shape)
    vault.add_damage_from_meshes(intrados, extrados)
    swt = vault.compute_selfweight()
    print('selfweight:', swt)
    if test_shape == 3:
        vault.data['b_manual'] = 0.0492494
    # view_shapes(vault).show()

    # ----------------------- 2. Create Form Diagram ---------------------------

    data_diagram = {
        'type': type_formdiagram,
        'H': H,
        'L': L,
        'total_nodes': discretisation,
        'x0': 0.0
    }

    form = FormDiagram.from_library(data_diagram)

    # --------------------- 3. Create Min Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'slsqp'
    optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
    optimiser.data['variables'] = ['ind', 'zb']
    optimiser.data['objective'] = 'min'  # Set the objective
    optimiser.data['printout'] = True
    optimiser.data['plot'] = False
    optimiser.data['find_inds'] = True
    optimiser.data['qmax'] = 1e+10  # Check if this is limiting the solution

    # --------------------- 4. Optimise ---------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope_with_damage()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    # ----------------------- 4. Save and retrieve parameters ---------------------------

    file_address = compas_tno.get('test.json')
    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight_min = fopt/swt
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight - min: ', fopt_over_weight_min)
    print('Exitflag: ', exitflag)
    if exitflag == 0:
        forms_min.append(deepcopy(form))
        shapes.append(vault)
        form.to_json(file_address)

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

    # ------------------------ 5 . Reduce the thickness ---------------------------

    if exitflag == 0:
        solutions_min.append(fopt_over_weight_min)
        solutions_max.append(fopt_over_weight_max)
        size_parameters.append(t_over_R)
        thk = round(thk - thk_reduction, 4)
        print('Reduce thickness to', thk, '\n')
    else:
        print('\nOptimisation did not find a solution for thk:', thk)
        print('Last solved thk:', round(thk + thk_reduction, 4))

# --------------- 6 . After exit print the list of solutions obtained ------------

print('\nSUMMARY')
print('\nSizes calculated')
print(size_parameters)
print('Solutions min')
print(solutions_min)
print('Solutions max')
print(solutions_max)

csv_path = os.path.join(compas_tno.get('/csv/'), 'diagram_{0}.csv'.format(test_shape))
save_csv(size_parameters, solutions_min, solutions_max, limit_state=True, path=csv_path)

img_graph = os.path.join(compas_tno.get('/imgs/'), 'diagram_dammage{0}.pdf'.format(test_shape))
diagram_of_thrust(size_parameters, solutions_min, solutions_max, xy_limits=[[0.20, 0.10], [60, 30]], fill=True, save=img_graph).show()
