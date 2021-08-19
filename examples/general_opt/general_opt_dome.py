from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution

from compas_tno.optimisers import Optimiser
from compas_tno.analysis.analysis import Analysis
import os

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

# Add this to the test

span = 10.0
k = 1.0
# discretisation = 10
type_formdiagram = 'radial_fd'
type_structure = 'dome'
thk = 0.50
radius = 5.0
span = radius * 2
n = 2

discretisation = [8, 12]

hc = None
he = None
lambd = 0.10

c = 0.05

save = False
solutions = {}

objective = ['min']
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'envelopexy']
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'
gradients = True

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])

print(variables)

for c in [c]:  # set the distance that the nodes can move
    solutions[c] = {}

    for obj in objective:  # set the objective
        solutions[c][obj] = {}

        for thk in [thk]:  # thickness of the problem

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'discretisation': discretisation,
                'center': [5.0, 5.0],
                # 'diagonal': True,
                # 'partial_diagonal': 'right',
                'radius': radius,
            }
            form = FormDiagram.from_library(data_diagram)

            # Create shape

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'discretisation': [discretisation[0]*n, discretisation[1]*n],
                'center': [5.0, 5.0],
                'radius': radius,
                't': 1.0,
            }
            vault = Shape.from_library(data_shape)
            vault.ro = 100.0

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            apply_envelope_from_shape(form, vault)
            apply_selfweight_from_shape(form, vault)
            if 'lambd' in variables:
                apply_horizontal_multiplier(form, lambd=lambd)

            if 'envelopexy' in constraints:
                apply_envelope_on_xy(form, c=c)
            apply_bounds_on_q(form, qmax=0.0)

            form_base = form.copy()

            # ------------------------------------------------------------
            # ------------------- Proper Implementation ------------------
            # ------------------------------------------------------------

            optimiser = Optimiser()
            optimiser.settings['library'] = solver
            optimiser.settings['solver'] = solver
            optimiser.settings['constraints'] = constraints
            optimiser.settings['variables'] = variables
            optimiser.settings['features'] = features
            optimiser.settings['axis_symmetry'] = axis_sym
            optimiser.settings['objective'] = obj
            optimiser.settings['plot'] = True
            optimiser.settings['find_inds'] = False
            optimiser.settings['printout'] = True
            optimiser.settings['max_iter'] = 2000
            optimiser.settings['gradient'] = gradients
            optimiser.settings['jacobian'] = gradients
            optimiser.settings['derivative_test'] = False

            optimiser.settings['starting_point'] = starting_point

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)
            # analysis.apply_selfweight()
            # analysis.apply_envelope()
            analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            analysis.run()

            form.overview_forces()
            if obj == 't':
                thk = form.attributes['thk']

            weight = 0
            for key in form.vertices():
                weight += form.vertex_attribute(key, 'pz')

            thrust = form.thrust()
            print('Ratio Thrust/Weight:', thrust/weight)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', type_structure, type_formdiagram)
            if 'ind' in optimiser.settings['variables']:
                folder = os.path.join(folder, 'fixed')
            else:
                folder = os.path.join(folder, 'mov_c_' + str(c))
            if he:
                folder = os.path.join(folder, 'hc_' + str(hc) + '_he_' + str(he))
            if 'lambd' in obj:
                folder = os.path.join(folder, 'max_lambd')
                lambd = -1 * optimiser.fopt
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)
            address = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'

            if 'lambd' in obj:
                address = save_form + '_' + optimiser.settings['objective'] + '_' + str(lambd) + '_thk_' + str(100*thk) + '.json'

            plot_superimposed_diagrams(form, form_base).show()
            plot_form(form, show_q=False, cracks=True).show()
            # view_solution(form).show()

            if optimiser.exitflag == 0:
                solutions[c][obj][thk] = thrust/weight * 100
                img_file = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.png'
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file).show()
                    plot_form(form, show_q=False, cracks=True).show()
            else:
                plot_superimposed_diagrams(form, form_base).show()
                view_solution(form).show()
                break

# form.to_json(address)
# print('Saved to: ', address)
from compas_tno.plotters.form import plot_form_semicirculararch_xz
from compas_tno.algorithms import reactions
reactions(form)
tol = 10e-3
form.attributes['Re'] = radius + thk/2
form.attributes['Ri'] = radius - thk/2
address_plot_section = os.path.join(folder, title) + '_' + obj + '_thk_' + str(thk) + '_plot_' + 'section' + '.pdf'
plot_form_semicirculararch_xz(form, radius=0.06, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True, yrange=[radius-tol, radius+tol], save=address_plot_section).show()

view_solution(form).show()

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])
