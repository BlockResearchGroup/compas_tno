import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.viewers.thrust import view_solution

# ----------------------------------------------------------------------------
# ---------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL FD --------------
# ----------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.50  # thickness on the start in meters
thk_reduction = 0.01  # in meters
solutions = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters

# Basic parameters

type_structure = 'dome'
type_formdiagram = 'spiral_fd'  # Try 'radial_spaced_fd' and 'spiral_fd'
discretisation = [8, 28]  # Try increasing a bit
R = 5.0

# ----------------------- 1. Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': R,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)
# plot_form(form).show()

# ----------------------- 3. Initiate loop on the optimisation ---------------------------

while exitflag == 0:

    t_over_R = thk/R
    print('\n----- Starting the' , type_structure, 'problem for thk/R:', t_over_R, '\n')

    # ------------------ 4. Create the new shape (extrados and intrados) for given thk ----------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': [2*discretisation[0], 2*discretisation[1]],
        'center': [5.0, 5.0],
        'radius': R,
        't' : 1.0
    }

    vault = Shape.from_library(data_shape)
    swt = vault.compute_selfweight()
    print('Selfweight computed:', swt)
    print('Dome created!')

    # --------------------- 5. Create Minimisation Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'slsqp'
    optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry']
    optimiser.data['variables'] = ['ind', 'zb']
    optimiser.data['objective'] = 'min'
    optimiser.data['printout'] = True
    optimiser.data['plot'] = False
    optimiser.data['find_inds'] = True
    optimiser.data['qmax'] = 3000.0
    print(optimiser.data)

    # --------------------- 6. Set up and run analysis ---------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    # --------------------- 7. Evaluate and plot results ---------------------

    file_address = compas_tno.get('test.json')
    form.to_json(file_address)
    exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
    fopt = optimiser.fopt  # objective function optimum value
    fopt_over_weight = fopt/swt  # divide by selfweight
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight: ', fopt_over_weight)
    print('Exitflag: ', exitflag)
    # plot_form(form, show_q=False, simple=True, cracks=True).show() # When cracks = True the intrados cracks will be blue and extrados green

    # ------------------------ 8 . Reduce the thickness ---------------------------

    if exitflag == 0:
        solutions.append(fopt_over_weight)
        size_parameters.append(t_over_R)
        thk = round(thk - thk_reduction, 4)
        print('Reduce thickness to', thk, '\n')
    else:
        print('\nOptimisation did not find a solution for thk:', thk)
        print('Last solved thk:', round(thk + thk_reduction, 4))


# --------------- 9 . After exit print the list of solutions obtained ------------

print('\n SUMMARY')
print('\ thk/R calculated')
print(size_parameters)
print('Solutions Found')
print(solutions)

view_thrust(form).show()

# If you wish to visualise the upper and lower bound together
# view_solution(form, vault).show()
