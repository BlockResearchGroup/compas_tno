from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis
import os
from compas_tno.plotters import save_csv
from compas_tno.plotters import diagram_of_thrust

span = 10.0
k = 1.0
# discretisation = 10
type_formdiagram = 'radial_fd'
type_structure = 'dome'
thk = 0.40
radius = 5.0
span = radius * 2
n = 2

discretisation = [20, 16]

hc = None
he = None
lambd = 0.10

c = 0.1

save = False
solutions = {}

objective = ['t']
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb']
features = ['fixed', 'sym']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])

print(variables)

for c in [0.1]:  # set the distance that the nodes can move
    solutions[c] = {}

    for obj in objective:  # set the objective
        solutions[c][obj] = {}

        for thk in [thk]:  # thickness of the problem

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'xy_span': [[0, span], [0, k*span]],
                'discretisation': discretisation,
                'center': [5.0, 5.0],
                'diagonal': False,
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

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            form.envelope_from_shape(vault)
            form.selfweight_from_shape(vault)
            if 'lambd' in variables:
                form.apply_horizontal_multiplier(lambd=lambd)

            # for key in form.vertices_where({'is_fixed': True}):
            #     form.vertex_attribute(key, 'z', 0.2)

            form.envelope_on_x_y(c=c)
            form.bounds_on_q(qmax=0.0)

            # address = '/Users/mricardo/compas_dev/me/general_opt/dome/radial_fd/mov_c_0.1/dome_radial_fd_discr_[20, 16]_t_thk_10.77604794596367.json'
            # form = FormDiagram.from_json(address)

            form_base = form.copy()

            # ------------------------------------------------------------
            # ------------------- Proper Implementation ------------------
            # ------------------------------------------------------------

            optimiser = Optimiser()
            optimiser.data['library'] = solver
            optimiser.data['solver'] = solver
            optimiser.data['constraints'] = constraints
            optimiser.data['variables'] = variables
            optimiser.data['features'] = features
            optimiser.data['axis_symmetry'] = axis_sym
            optimiser.data['objective'] = obj
            optimiser.data['plot'] = False
            optimiser.data['find_inds'] = False
            optimiser.data['printout'] = True
            optimiser.data['max_iter'] = 2000
            optimiser.data['gradient'] = True
            optimiser.data['jacobian'] = True
            optimiser.data['derivative_test'] = True

            optimiser.data['starting_point'] = starting_point

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
            if 'ind' in optimiser.data['variables']:
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
            address = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'

            if 'lambd' in obj:
                address = save_form + '_' + optimiser.data['objective'] + '_' + str(lambd) + '_thk_' + str(100*thk) + '.json'

            plot_superimposed_diagrams(form, form_base).show()
            plot_form(form, show_q=False, cracks=True).show()
            # view_solution(form).show()

            if optimiser.exitflag == 0:
                solutions[c][obj][thk] = thrust/weight * 100
                img_file = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.png'
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file).show()
                    plot_form(form, show_q=False, cracks=True).show()
            else:
                plot_superimposed_diagrams(form, form_base).show()
                view_solution(form).show()
                break

    form.to_json(address)
    print('Saved to: ', address)
    from compas_tno.plotters.form import plot_form_semicirculararch_xz
    from compas_tno.algorithms import reactions
    reactions(form, plot=True)
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
