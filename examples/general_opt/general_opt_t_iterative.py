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
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk0 = 0.50
discretisation_shape = 10 * discretisation


hc = None
he = None
plot_comparison = True
update_loads = True
lambd = None

c = 0.1

thk_reduction = 0.05
save = False
solutions = {}

solver = 'IPOPT'
constraints = ['envelope', 'envelopexy']
variables = ['q', 'zb', 't']
features = ['sym', 'fixed']
starting_point = 'loadpath'

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
            # vault.ro = 20.0

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
            # apply_bounds_on_q(form, qmax=0.0)
            apply_bounds_on_q(form)

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
            optimiser.settings['objective'] = obj
            optimiser.settings['plot'] = False
            optimiser.settings['find_inds'] = False
            optimiser.settings['max_iter'] = 1000
            optimiser.settings['gradient'] = True
            optimiser.settings['printout'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['derivative_test'] = False
            optimiser.settings['starting_point'] = starting_point

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)

            min_thks = []
            errors = []
            entropies = []
            entropiespz = []

            entropy = 1.0
            entropypz = 0.0

            pz0 = array(form.vertices_attribute('pz')).reshape(-1, 1)

            i = 0

            while entropy > 10e-4 or i < 2:

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

                solutions[c][obj][thk] = thrust/weight * 100
                img_file = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.png'
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file)

                vault.data['thk'] = thk

                apply_envelope_from_shape(form, vault)
                if update_loads:  # Update bounds and SWT
                    apply_selfweight_from_shape(form, vault)

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

                plot_superimposed_diagrams(form, form_base).show()
                # plot_form(form, show_q=False, cracks=True).show()
                # view_solution(form).show()

                optimiser.settings['starting_point'] = 'current'
                optimiser.settings['max_thk'] = min(thk*1.25, 0.5)
                apply_bounds_on_q(form, qmax=0.0)
                optimiser.settings['features'] = ['sym', 'adapted-envelope']
                i += 1

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
