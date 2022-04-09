from compas_tno.diagrams import FormDiagram
from compas_tno.plotters.plotter import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from compas_view2.shapes import Arrow

from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import norm_vector
from numpy import array

import os
import math

span = 10.0
k = 1.0
discretisation = 14
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 2 * discretisation
hc = None
he = None

corner = True

save = True
solutions = {}

objective = ['Ecomp-linear']
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']  # ['sym', 'update-envelope']
axis_sym = [[[0.0, 0.0], [10.0, 10.0]]]
# qmax = 10e+6
starting_point = 'loadpath'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambdh']:
    variables.append(objective[0])
    lambd = 0.1

# min_thk_json = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_14_min_thk_t_0.3355366629102956.json'

Xc = [5.0, 5.0, 0.0]
c = 0.1
phi = 0.0

phi_list = [0.0]  # [0.0]  # [0.0, 22.5, 45, 67.5, 90.0]  # [-22.5, -45, -67.5]
ro = 45.0

for phi in phi_list:
    solutions[phi] = {}
    for sign in [1]:  # get the minimum and maximum
        solutions[phi][sign] = {}

        for obj in objective:  # set the objective
            solutions[phi][sign][obj] = {}

            for thk in [0.50]:  # thickness of the problem # [0.50, 0.45, 0.40, 0.35, 0.336]

                # min_thk_json = '/Users/mricardo/compas_dev/me/min_thk/crossvault/cross_fd/crossvault_cross_fd_discr_14_min_thk_t_0.3355366629102956.json'

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
                vault.ro = 20.0

                solutions[phi][sign][obj][thk] = {}

                # ------------------------------------------------------------
                # -----------------------  INITIALISE   ----------------------
                # ------------------------------------------------------------

                # Apply Selfweight and Envelope

                apply_envelope_from_shape(form, vault)
                apply_selfweight_from_shape(form, vault)

                if 'envelopexy' in constraints:
                    apply_envelope_on_xy(form, c=c)
                apply_bounds_on_q(form)

                form_base = form.copy()

                # Consider the displacement vecto

                lines = []
                vector_supports = []

                # sign = +1  # +1 for outwards / -1 for inwards

                view = Viewer(form, vault)
                vectors = []
                bases = []

                if corner:

                    for key in form.vertices_where({'is_fixed': True}):
                        x, y, z = form.vertex_coordinates(key)

                        if x > Xc[0] and y > Xc[1]:             # 1 corner only
                            # dXbi = scale_vector(normalize_vector([1, 1, math.sin(math.rad(phi))]), sign)

                            dXbi = [sign * math.cos(math.radians(ro)) * math.cos(math.radians(-phi)),
                                    sign * math.cos(math.radians(ro)) * math.cos(math.radians(-phi)),
                                    sign * math.sin(math.radians(-phi))]
                            print('Norm of vector:', norm_vector(dXbi))
                        else:
                            dXbi = [0, 0, 0]

                        vector_supports.append(dXbi)
                        lines.append({
                            'start': [x, y, z],
                            'end': [x + dXbi[0], y + dXbi[1], z + dXbi[2]],
                            'width': 3
                        })

                        if norm_vector(dXbi) > 0.01:
                            arrow = Arrow([x, y, z], dXbi)
                            bases.append([x, y, z])
                            vectors.append(arrow)
                            view.app.add(arrow, color=(0, 0, 0))

                    dXb = array(vector_supports)
                    print(dXb)

                    view.draw_shape()
                    view.show()

                else:
                    raise NotImplementedError

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
                optimiser.settings['max_iter'] = 500
                optimiser.settings['gradient'] = True
                optimiser.settings['jacobian'] = True
                optimiser.settings['printout'] = True
                optimiser.settings['jacobian'] = True
                optimiser.settings['derivative_test'] = False
                optimiser.settings['starting_point'] = starting_point
                optimiser.settings['axis_symmetry'] = axis_sym
                optimiser.settings['support_displacement'] = dXb
                optimiser.settings['Ecomp_method'] = 'complete'

                # --------------------- 5. Set up and run analysis ---------------------

                analysis = Analysis.from_elements(vault, form, optimiser)
                analysis.set_up_optimiser()
                analysis.run()

                form.overview_forces()
                if obj == 't':
                    thk = form.attributes['thk']
                if 'Ecomp' in obj.split('-'):
                    Ecomp = sign*optimiser.fopt

                weight = 0
                for key in form.vertices():
                    weight += form.vertex_attribute(key, 'pz')

                thrust = form.thrust()

                T_over_W = abs(thrust/weight)
                E_over_W = abs(Ecomp/weight)

                print('Ratio Thrust/Weight:', T_over_W)
                print('Ratio Ecomp/Weight:', E_over_W)

                folder = os.path.join('/Users/mricardo/compas_dev/me', 'compl_energy', type_structure, type_formdiagram)
                if 'fixed' in optimiser.settings['features']:
                    folder = os.path.join(folder, 'fixed')
                else:
                    folder = os.path.join(folder, 'mov_c_' + str(c))
                    if 'update-envelope' in optimiser.settings['features']:
                        folder = os.path.join(folder, 'update-envelope')
                if he:
                    folder = os.path.join(folder, 'hc_' + str(hc) + '_he_' + str(he))
                if corner:
                    folder = os.path.join(folder, 'corner')
                if phi:
                    folder = os.path.join(folder, 'phi_' + str(phi))
                if sign:
                    folder = os.path.join(folder, 'sign_' + str(sign))
                os.makedirs(folder, exist_ok=True)
                title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
                save_form = os.path.join(folder, title)
                address = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'

                # plot_superimposed_diagrams(form, form_base).show()
                # view = Viewer(form)
                # view.show_solution()

                print('Optimiser exitflag:', optimiser.exitflag)

                if optimiser.exitflag == 0:
                    img_file = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.png'
                    solutions[phi][sign][obj][thk]['T/W'] = T_over_W
                    solutions[phi][sign][obj][thk]['E/W'] = E_over_W
                    if save:
                        form.to_json(address)
                        print('Saved to: ', address)
                        # plot_superimposed_diagrams(form, form_base, save=img_file)#.show()
                        # plot_form(form, show_q=False, cracks=True, save=img_file)  #.show()
                        # starting_point = 'current'
                else:
                    # plot_superimposed_diagrams(form, form_base).show()
                    # view = Viewer(form)
                    # view.show_solution()
                    solutions[phi][sign][obj][thk]['T/W'] = 'ERROR'
                    solutions[phi][sign][obj][thk]['E/W'] = 'ERROR'
                    starting_point = 'loadpath'

                print(solutions)
                print('\n')

view = Viewer(form)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
for vector in vectors:
    view.app.add(vector, color=(0, 0, 0))
view.show()

plotter = TNOPlotter(form, shape=vault)
plotter.draw_form()
plotter.draw_cracks()
plotter.draw_vectors(vectors=[vector.direction for vector in vectors], bases=[vector.position for vector in vectors])
# plotter.draw_reactions()
plotter.show()

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])
