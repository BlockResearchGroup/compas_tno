import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis.analysis import Analysis

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR ARCH --------------------
# ----------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.200  # in meters
thk_reduction = 0.005  # in meters
solutions = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters

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
    print('\n----- The ', type_formdiagram, 'problem for thk/R:', t_over_R, '\n')

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
    fopt = optimiser.fopt
    fopt_over_weight = fopt/swt
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight: ', fopt_over_weight)
    print('Exitflag: ', exitflag)
    plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()

    # ------------------------ 5 . Reduce the thickness ---------------------------

    if exitflag == 0:
        solutions.append(fopt_over_weight)
        size_parameters.append(t_over_R)
        thk = round(thk - thk_reduction, 4)
        print('Reduce thickness to', thk, '\n')
    else:
        print('\nOptimisation did not find a solution for thk:', thk)
        print('Last solved thk:', round(thk + thk_reduction, 4))

# --------------- 6 . After exit print the list of solutions obtained ------------

print('\n SUMARY')
print('\nSizes calculated')
print(size_parameters)
print('Solutions Found')
print(solutions)
