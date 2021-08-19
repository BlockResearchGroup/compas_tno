import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis.analysis import Analysis
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_multiple_thrust
from compas_tno.plotters import surface_GSF_load_mult

# ----------------------------------------------------------------------
# ---- EXAMPLE OF MIN and MAX THRUST FOR ARCH WITH INCREMENTAL THK -----
# ----------------------------------------------------------------------

# Setup Parameters

exitflag = 0  # means that optimisation found a solution
thk = 0.200  # thickness on the start in meters
thk_reduction = 0.001  # in meters
solutions_min = []  # empty lists to keep track of  the solutions
solutions_max = []  # empty lists to keep track of  the solutions
size_parameters = []  # empty lists to keep track of  the parameters
hor_mult = 0.10
direction_loads = 'px'


# Basic parameters

H = 1.0
L = 2.0
R = H/2 + (L**2 / (8 * H))  # Note that this reduces to L/2 when springing angle = 180 deg, or L = 2*H
discretisation = 20
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


# --------------------- 2. Create Optimiser still with no objective  ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 10000.0  # Check if this is limiting the solution

hor_study = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.293]
sizes = []
maxs = []
mins = []
legends = []
legends_num = []
for hor_mult in hor_study:

    exitflag = 0  # means that optimisation found a solution
    thk = 0.200  # thickness on the start in meters
    solutions_min = []  # empty lists to keep track of  the solutions
    solutions_max = []  # empty lists to keep track of  the solutions
    size_parameters = []  # empty lists to keep track of  the parameters
    while exitflag == 0:

        t_over_R = thk/R
        print('\n----- The' , type_formdiagram, 'problem for thk/R:', t_over_R, '\n')

        # -------------------------- 1. Create Arch shape for thk ---------------------------

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
        analysis.apply_hor_multiplier(hor_mult, direction_loads)

        # --------------------------- 4.1 Set The objective to min and run ---------------------------

        optimiser.settings['objective'] = 'min'
        analysis.set_up_optimiser()
        analysis.run()

        # --------------------------- 4.2 Extract Data from Optimisation ---------------------------

        exitflag = optimiser.exitflag  # get info if optimisation was succeded ot not
        fopt = optimiser.fopt  # objective function optimum value
        fopt_over_weight_min = fopt/swt
        print('Thickness', thk)
        print('thk/R:', t_over_R)
        print('Thrust over weight - min: ', fopt_over_weight_min)
        print('Exitflag: ', exitflag)
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
        print('thk/R:', t_over_R)
        print('Thrust over weight - max: ', fopt_over_weight_max)
        print('Exitflag: ', exitflag)
        # plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()

        # --------------------------- 6. Reduce the thickness ---------------------------

        if exitflag == 0:
            solutions_min.append(fopt_over_weight_min)
            solutions_max.append(fopt_over_weight_max)
            size_parameters.append(t_over_R)
            thk = round(thk - thk_reduction, 4)
            print('Reduce thickness to', thk, '\n')
        else:
            print('\nOptimisation did not find a solution for thk:', thk)
            print('Last solved thk:', round(thk + thk_reduction, 4))

    if size_parameters != []:
        # diagram_of_thrust(size_parameters, solutions_min, solutions_max).show()
        maxs.append(solutions_max)
        mins.append(solutions_min)
        sizes.append(size_parameters)
        legends.append(str('px='+str(hor_mult)))
        legends_num.append(hor_mult)


# --------------- 6 . After exit print the list of solutions obtained ------------

print('\nSUMMARY')

print('sizes')
print(sizes)
print('mins')
print(mins)
print('maxs')
print(maxs)
print('legend')
print(legends)
# diagram_of_thrust(size_parameters, solutions_min, solutions_max).show()
# diagram_of_multiple_thrust(sizes, mins, maxs, legends, limit_state=False).show()
surface_GSF_load_mult(sizes, mins, maxs, legends_num).show()

