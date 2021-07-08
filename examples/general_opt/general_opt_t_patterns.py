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
save = True
solutions = {}
hc = 6.0
he = 5.0

# for c in [0.1, 0.25, 0.50]:
for c in [0.1]:
    solutions[c] = {}

    for obj in ['t']:
        solutions[c][obj] = {}

        # for thk in [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]:
        # for type_ in ['B2', 'D3', 'D2', 'D4']:
        # for type_ in ['B1', 'B2', 'C1', 'C2', 'D2', 'D3', 'D4']:
        for type_ in ['B3', 'C1', 'C2', 'D2']:

            solutions[c][obj][type_] = {}

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

            # Create form diagram

            thk = 0.50

            pattern_fd = '/Users/mricardo/compas_dev/me/loadpath/corner/topology/' + type_ + '.json'

            form = FormDiagram.from_json(pattern_fd)
            # plot_form(form, show_q=False).show()

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            form.envelope_from_shape(vault)
            form.selfweight_from_shape(vault)

            form.envelope_on_x_y(c=c)

            form_base = form.copy()

            # apply_sag(form)

            # folder_lp = os.path.join('/Users/mricardo/compas_dev/me', 'loadpath', type_structure, type_formdiagram)
            # os.makedirs(folder_lp, exist_ok=True)
            # title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            # save_lp = os.path.join(folder_lp, title)
            # address_lp = save_lp + '_' + 'lp' + '_thk_' + str(100*thk) + '.json'
            # try:
            #     form = FormDiagram.from_json(save_lp)
            #     form.envelope_from_shape(vault)
            #     form.selfweight_from_shape(vault)
            #     form.envelope_on_x_y(c=c)
            # except:
            #     form.initialise_loadpath()
            #     form.to_json(address_lp)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', type_structure, type_, 'mov_c_' + str(c), 'hc_' + str(hc) + '_he_' + str(he))
            os.makedirs(folder, exist_ok=True)
            # title = type_structure + '_' + type_ + '_discr_' + str(discretisation)
            # save_form = os.path.join(folder, title)
            # address_load = save_form + '_' + 'min' + '_thk_' + str(100*thk) + '.json'
            # form = FormDiagram.from_json(address_load)

            # form.envelope_on_x_y(c=c)

            # plot_form(form, show_q=False, cracks=True).show()
            # view_solution(form).show()

            # ------------------------------------------------------------
            # ------------------- Proper Implementation ------------------
            # ------------------------------------------------------------

            optimiser = Optimiser()
            # optimiser.data['library'] = 'SLSQP'
            # optimiser.data['solver'] = 'SLSQP'
            optimiser.data['library'] = 'IPOPT'
            optimiser.data['solver'] = 'IPOPT'
            optimiser.data['constraints'] = ['funicular', 'envelope']
            optimiser.data['variables'] = ['sym', 'zb', 't']
            # optimiser.data['variables'] = ['ind', 'zb', 't']
            optimiser.data['objective'] = obj
            optimiser.data['plot'] = False
            optimiser.data['find_inds'] = False
            optimiser.data['max_iter'] = 500
            optimiser.data['qmax'] = 1000.0
            optimiser.data['gradient'] = True
            optimiser.data['printout'] = True
            optimiser.data['jacobian'] = True
            optimiser.data['derivative_test'] = False
            optimiser.data['max_iter'] = 500

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)
            # analysis.apply_selfweight()
            # analysis.apply_envelope()
            # analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            analysis.run()

            form.overview_forces()
            thk = form.attributes['thk']

            weight = 0
            for key in form.vertices():
                weight += form.vertex_attribute(key, 'pz')

            thrust = form.thrust()
            print('Ratio Thrust/Weight:', thrust/weight)

            # folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', 'min_thk', type_structure, type_formdiagram, 'mov_c_' + str(c))
            # os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_
            save_form = os.path.join(folder, title)
            address = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'
            img_file = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.png'

            # plot_superimposed_diagrams(form, form_base).show()
            # view_solution(form).show()

            if optimiser.exitflag == 0:
                solutions[c][obj][thk] = thrust/weight * 100
                solutions[c][obj][type_] = thk
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file)
                    # view_solution(form).show()
            else:
                plot_superimposed_diagrams(form, form_base, save=img_file)
                # break

    # plot_superimposed_diagrams(form, form_base).show()
    # view_solution(form).show()

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
