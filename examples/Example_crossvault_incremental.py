import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.viewers.thrust import view_solution

# ----------------------------------------------------------------------------
# ------ EXAMPLE OF INCREMENTAL MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ----------------------------------------------------------------------------


exitflag = 0  # means that optimisation found a solution
thk = 0.50  # thickness on the start in meters
thk_reduction = 0.01  # in meters
solutions = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
span = 10.0  # square span for analysis
R = span/2  # Only valid for rounded cross vault

# Basic parameters

type_structure = 'crossvault'
type_formdiagram = 'cross_fd'  # Try also 'fan_fd'
discretisation = 10

# ----------------------- 1. Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 2. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)
plot_form(form).show()

# ----------------------- 3. Initiate loop on the optimisation ---------------------------

while exitflag == 0:

    t_over_R = thk/R
    print('\n----- Starting the' , type_formdiagram, 'problem for thk/R:', t_over_R, '\n')

    # ------------------ 4. Create the new shape (extrados and intrados) for given thk ----------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation,
        'xy_span': [[0, span], [0, span]],
        't': 0.0
    }

    vault = Shape.from_library(data_shape)
    swt = vault.compute_selfweight()
    print('Selfweight computed:', swt)
    print('Crossvault created!')

    # --------------------- 5. Create Minimisation Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'slsqp'
    optimiser.data['constraints'] = ['funicular', 'envelope']
    optimiser.data['variables'] = ['ind', 'zb']
    optimiser.data['objective'] = 'min'
    optimiser.data['printout'] = True
    optimiser.data['plot'] = False
    optimiser.data['find_inds'] = True
    optimiser.data['qmax'] = 2000.0
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
    fopt_over_weight = fopt/swt # divide by selfweight
    print('Thickness', thk)
    print('thk/R:', t_over_R)
    print('Thrust over weight: ', fopt_over_weight)
    print('Exitflag: ', exitflag)
    plot_form(form, show_q=False, simple=True, cracks=True).show() # When cracks = True the intrados cracks will be blue and extrados green

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
