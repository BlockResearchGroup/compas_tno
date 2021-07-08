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
discretisation = 14
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 10 * discretisation
hc = None
he = None

c = 0.1

save = False
solutions = {}

objective = ['t']
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 't']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])

for c in [0.1]:  # set the distance that the nodes can move
    solutions[c] = {}

    for obj in objective:  # set the objective
        solutions[c][obj] = {}

        for thk in [0.50]:  # thickness of the problem

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'xy_span': [[0, span], [0, k*span]],
                'discretisation': discretisation,
                'fix': 'corners'
            }

            form = FormDiagram.from_library(data_diagram)

            # Create shape

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'discretisation': discretisation_shape,
                'xy_span': [[0, span], [0, k*span]],
                'hc': hc,
                'hm': None,
                'he': he if he is None else [he, he, he, he],
                'center': [5.0, 5.0],
                'radius': span/2,
                't': 0.0,
            }

            vault = Shape.from_library(data_shape)

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            form.envelope_from_shape(vault)
            form.selfweight_from_shape(vault)

            form.envelope_on_x_y(c=c)
            form.bounds_on_q()

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
            optimiser.data['objective'] = obj
            optimiser.data['plot'] = False
            optimiser.data['find_inds'] = False
            optimiser.data['max_iter'] = 500
            optimiser.data['gradient'] = True
            optimiser.data['jacobian'] = True
            # optimiser.data['starting_point'] = 'loadpath'
            optimiser.data['printout'] = True
            optimiser.data['jacobian'] = True
            optimiser.data['derivative_test'] = True

            optimiser.data['starting_point'] = starting_point

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)
            # analysis.apply_selfweight()
            # analysis.apply_envelope()
            # analysis.apply_reaction_bounds()
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

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', 'min_thk', type_structure, type_formdiagram)
            if 'ind' in optimiser.data['variables']:
                folder = os.path.join(folder, 'fixed')
            else:
                folder = os.path.join(folder, 'mov_c_' + str(c))
            if he:
                folder = os.path.join(folder, 'hc_' + str(hc) + '_he_' + str(he))
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)
            address = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'

            plot_superimposed_diagrams(form, form_base).show()
            view_solution(form).show()

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

    view_solution(form).show()


print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])

# ----------------------- 5. Create Analysis loop on limit analysis --------------------------

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', type_structure, type_formdiagram, 'mov_c_' + str(c))
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# save_form = os.path.join(folder, title)

# analysis = Analysis.from_elements(vault, form, optimiser)
# results = analysis.limit_analysis_GSF(thk, thk_reduction, span, save_forms=save_form)
# thicknesses, size_parameters, solutions_min, solutions_max = results

# # ----------------------- 6. Save output data --------------------------

# csv_file = os.path.join(folder, title + '_data.csv')
# save_csv(size_parameters, solutions_min, solutions_max, path=csv_file, title=title)

# img_graph = os.path.join(folder, title + '_diagram.pdf')
# diagram_of_thrust(size_parameters, solutions_min, solutions_max, save=img_graph).show()
