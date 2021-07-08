from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution

import time

from compas_tno.algorithms import apply_sag

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis
import os
from compas_tno.plotters import save_csv
from compas_tno.plotters import diagram_of_thrust

from numpy import array

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'fan_fd'
type_structure = 'crossvault'
thk0 = 0.50
discretisation_shape = 10 * discretisation


hc = None
he = None
plot_comparison = True
update_loads = False

c = 0.1

thk_reduction = 0.05
save = False
solutions = {}

solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'envelopexy']
variables = ['q', 'zb', 't']
features = ['sym']

t0 = time.time()

# for c in [0.1, 0.25, 0.50]:
for c in [0.1]:
    solutions[c] = {}

    for obj in ['t']:
        solutions[c][obj] = {}

        # for thk in [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]:
        for thk in [thk0]:

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
                # 'hc': hc,
                # 'hm': None,
                # 'he': [he, he, he, he],
                'center': [5.0, 5.0],
                'radius': span/2,
                't': 0.0,
            }

            vault = Shape.from_library(data_shape)

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            form.envelope_from_shape(vault)
            form.selfweight_from_shape(vault)
            form.bounds_on_q()

            pz0 = array(form.vertices_attribute('pz')).reshape(-1, 1)

            form.envelope_on_x_y(c=c)

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
            optimiser.data['max_iter'] = 1500
            optimiser.data['gradient'] = True
            optimiser.data['printout'] = True
            optimiser.data['jacobian'] = True
            optimiser.data['derivative_test'] = False

            optimiser.data['starting_point'] = 'sag'

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)

            min_thks = []
            errors = []
            entropies = []
            entropiespz = []

            entropy = 1.0
            entropypz = 0.0

            while entropy > 10e-4:

                analysis.set_up_optimiser()
                analysis.run()

                form.overview_forces()
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

                solutions[c][obj][thk] = thrust/weight * 100
                img_file = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.png'
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file)

                vault.data['thk'] = thk

                form.envelope_from_shape(vault)
                if update_loads:  # Update bounds and SWT
                    form.selfweight_from_shape(vault)

                x = array(form.vertices_attribute('x')).reshape(-1, 1)
                y = array(form.vertices_attribute('y')).reshape(-1, 1)
                pzf = array(form.vertices_attribute('pz')).reshape(-1, 1)

                diff_xy = (x - optimiser.M.x0) ** 2 + (y - optimiser.M.y0) ** 2
                entropy = float(sum(diff_xy))
                entropypz = sum((pz0 - pzf)**2)

                optimiser.M.x0 = array(form.vertices_attribute('x')).reshape(-1, 1)
                optimiser.M.y0 = array(form.vertices_attribute('y')).reshape(-1, 1)
                optimiser.x0 = optimiser.xopt
                gfinal = optimiser.fconstr(optimiser.x0, optimiser.M)
                print('Max / Min g final:', max(gfinal), min(gfinal))
                print('Entropy XY:', entropy)
                print('Entropy Pz:', entropypz)

                pz0 = pzf

                entropies.append(entropy)
                entropiespz.append(entropypz)
                constr_error = 0.0
                for el in gfinal:
                    if el < 0:
                        constr_error += el**2
                print('Min thk:', thk)
                print('Error observed:', constr_error)
                min_thks.append(thk)
                errors.append(constr_error)
                # plot_form(form, show_q=False, cracks=True).show()

                # plot_superimposed_diagrams(form, form_base).show()
                # plot_form(form, show_q=False, cracks=True).show()
                # view_solution(form).show()

                optimiser.data['starting_point'] = 'current'
                optimiser.data['max_thk'] = min(thk*1.25, 0.5)

print(solutions)
print('\n')

print(min_thks)
print(entropies)
print(errors)

if plot_comparison:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    xaxis1 = list(range(len(min_thks) + 1))
    xaxis2 = xaxis1[1:]

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('iteration')
    ax1.set_ylabel('minimum thickness [m]')
    ln1 = ax1.plot(xaxis1, [thk0] + min_thks, label='minimum thickness', color='C1')

    ax2 = ax1.twinx()
    ax2.set_ylabel(r'horizontal nodal movement $[\mathrm{m}^2]$')
    ln2 = ax2.plot(xaxis2, entropies, label='entropy', color='C2')

    lns = ln1+ln2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)

    fig.tight_layout()
    plt.show()

plot_superimposed_diagrams(form, form_base).show()

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])

elapsed_time = time.time() - t0
print('Tipe elapsed in Script:', time)

view_solution(form).show()

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
