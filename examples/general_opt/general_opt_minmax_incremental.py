from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution

from compas_tno.algorithms import apply_sag

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis
import os
from compas_tno.plotters import save_csv
from compas_tno.plotters import diagram_of_thrust

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
type_structure = 'pointed_crossvault'
thk = 0.50
discretisation_shape = 4 * discretisation

c = 0.1

thk = 0.50
thk_reduction = 0.05

solutions = {}
sag = False
start_from_form = True
save = True

hc = 6.0
he = 5.0

adress = None

# address_solved = '/Users/mricardo/compas_dev/me/general_opt/min_thk/crossvault/fan_fd/mov_c_0.5/crossvault_fan_fd_discr_10_t_thk_17.63944480658616.json'
# address_solved = '/Users/mricardo/compas_dev/me/general_opt/min_thk/crossvault/fan_fd/mov_c_0.5/crossvault_fan_fd_discr_10_t_thk_13.340323711187377.json'
# address_solved = '/Users/mricardo/compas_dev/me/general_opt/crossvault/fan_fd/mov_c_0.5/crossvault_fan_fd_discr_10_min_thk_15.0.json'

for c in [0.1]:
    solutions[c] = {}

    for obj in ['min', 'max']:
        solutions[c][obj] = {}

        for thk in [0.50]:

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'xy_span': [[0, span], [0, k*span]],
                'discretisation': discretisation,
                'fix': 'corners'
            }

            form = FormDiagram.from_library(data_diagram)
            # plot_form(form, show_q=False).show()

            # Create shape

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'discretisation': discretisation_shape,
                'xy_span': [[0, span], [0, k*span]],
                'hc': hc,
                'hm': None,
                'he': [he, he, he, he],
                'center': [5.0, 5.0],
                'radius': span/2,
                't': 0.0,
            }

            vault = Shape.from_library(data_shape)

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            from compas_tno.utilities import apply_envelope_from_shape
            from compas_tno.utilities import apply_selfweight_from_shape
            from compas_tno.utilities import apply_envelope_on_xy
            from compas_tno.utilities import apply_horizontal_multiplier
            from compas_tno.utilities import apply_bounds_on_q

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
            # optimiser.settings['library'] = 'SLSQP'
            # optimiser.settings['solver'] = 'SLSQP'
            optimiser.settings['library'] = 'IPOPT'
            optimiser.settings['solver'] = 'IPOPT'
            optimiser.settings['constraints'] = ['funicular', 'envelope']
            optimiser.settings['variables'] = ['q', 'zb']
            # optimiser.settings['variables'] = ['sym', 'zb']
            # optimiser.settings['variables'] = ['q', 'zb']
            optimiser.settings['objective'] = obj
            optimiser.settings['plot'] = False
            optimiser.settings['find_inds'] = False
            optimiser.settings['qmax'] = 1000.0
            optimiser.settings['gradient'] = True
            optimiser.settings['printout'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['derivative_test'] = False
            optimiser.settings['max_iter'] = 1000

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)
            # analysis.apply_selfweight()
            # analysis.apply_envelope()
            # analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            analysis.run()

            weight = 0
            for key in form.vertices():
                weight += form.vertex_attribute(key, 'pz')

            form.overview_forces()
            thrust = abs(optimiser.fopt)
            print('Ratio Thrust/Weight:', thrust/weight)

            # plot_form(form).show()
            # view_solution(form).show()

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', 'min_thk', type_structure, type_formdiagram)
            if 'ind' in optimiser.settings['variables']:
                folder = os.path.join(folder, 'fixed')
            else:
                folder = os.path.join(folder, 'mov_c_' + str(c))
            if he:
                folder = os.path.join(folder, 'hc_' + str(hc) + '_he_' + str(he))
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)
            address = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'

            if optimiser.exitflag == 0:
                solutions[c][obj][thk] = thrust/weight * 100
                if save:
                    form.to_json(address)
                    print('Form saved to:', address)
                    img_file = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.png'
                    plot_superimposed_diagrams(form, form_base, save=img_file)
                else:
                    plot_superimposed_diagrams(form, form_base).show()
                    view_solution(form).show()
            else:
                plot_superimposed_diagrams(form, form_base)
                pass

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        for key3 in solutions[key][key2]:
            print(key, key2, key3, solutions[key][key2][key3])

# # # ----------------------- 5. Create Analysis loop on limit analysis --------------------------

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
