from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from compas.geometry import Vector
from compas.geometry import Point

from compas.geometry import normalize_vector
from numpy import array
import os

span = 10.0
k = 1.0
discretisation = [20, 16]
type_formdiagram = 'radial_fd'
type_structure = 'dome'
thk = 0.50
discretisation_shape = [2 * discretisation[0], 2 * discretisation[1]]

save = True
solutions = {}

objective = ['Ecomp-linear']
Ecomp_method = 'complete'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']  # , 'envelopexy'
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None
# axis_sym = [[0.0, 5.0], [10.0, 5.0]]
# axis_sym = [[5.0, 0.0], [5.0, 10.0]]
axis_sym = None
# qmax = 10e+6
starting_point = 'loadpath'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])
    lambd = 0.1

Xc = [5.0, 5.0, 0.0]

for c in [0.1]:  # set the distance that the nodes can move
    solutions[c] = {}

    for obj in objective:  # set the objective
        solutions[c][obj] = {}

        for thk in [0.50]:  # thickness of the problem

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'center': [5.0, 5.0],
                'radius': span/2,
                'discretisation': discretisation
            }

            form = FormDiagram.from_library(data_diagram)

            # Create shape

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'discretisation': discretisation_shape,
                'center': [5.0, 5.0],
                'radius': span/2,
                't': 1.0,
            }

            vault = Shape.from_library(data_shape)
            vault.ro = 20.0

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

            # Consider the displacement vector

            if 'Ecomp' in obj.split('-'):

                lines = []
                vector_supports = []
                vectors = []
                bases = []

                sign = -1  # +1 for outwards / -1 for inwards

                for key in form.vertices_where({'is_fixed': True}):
                    x, y, z = form.vertex_coordinates(key)

                    # dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), sign*(z - Xc[2])])

                    if abs(Xc[0] - x) > 10e-3:
                        dXbi = normalize_vector([sign*(Xc[0] - x), 0, 0])  # vault opening up
                    else:
                        dXbi = [0, 0, 0]

                    # dXbi = normalize_vector([1, 0, 0])  # rigid body

                    # if x > Xc[0]:         # If else to move 2 supports left
                    #     dXbi = normalize_vector([1, 0, 0])
                    # else:
                    #     dXbi = [0, 0, 0]

                    # if x > Xc[0] and y > Xc[1]:             # 1 corner only
                    #     dXbi = normalize_vector([1, 1, 0])
                    # else:
                    #     dXbi = [0, 0, 0]

                    # if x > Xc[0] and y > Xc[1]:             # vertical settlement
                    #     dXbi = normalize_vector([0, 0, -1])
                    # else:
                    #     dXbi = [0, 0, 0]

                    # if x - Xc[0] > 2.5:
                    #     dXbi = normalize_vector([0, 0, -1])  # vertical settlement right end
                    # else:
                    #     dXbi = [0, 0, 0]

                    vector_supports.append(dXbi)
                    vectors.append(Vector(*dXbi))
                    bases.append(Point(x, y, z))
                    # lines.append({
                    #     'start': [x, y, z],
                    #     'end': [x + dXbi[0], y + dXbi[1], z + dXbi[2]],
                    #     'width': 3
                    # })

                dXb = array(vector_supports)
                print(dXb)

                from compas_tno.plotters import TNOPlotter

                plotter = TNOPlotter(form)
                plotter.draw_mesh()
                plotter.draw_vectors(vectors, bases)
                plotter.show()

                # key_index = form.key_index()
                # plotter = MeshPlotter(form)
                # plotter.draw_edges()
                # plotter.draw_vertices(keys=form.fixed(), facecolor={key: '000000' for key in form.vertices_where({'is_fixed': True})})
                # plotter.draw_arrows(lines)
                # plotter.show()

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
            optimiser.settings['axis_symmetry'] = axis_sym
            optimiser.settings['max_iter'] = 500
            optimiser.settings['gradient'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['printout'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['derivative_test'] = False
            optimiser.settings['starting_point'] = starting_point
            optimiser.settings['support_displacement'] = dXb
            optimiser.settings['Ecomp_method'] = Ecomp_method
            optimiser.settings['save_iterations'] = True

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

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'compl_energy', 'assym', type_structure, type_formdiagram)
            if 'ind' in optimiser.settings['variables']:
                folder = os.path.join(folder, 'fixed')
            else:
                if 'fixed' not in features:
                    folder = os.path.join(folder, 'mov_c_' + str(c))
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
                solutions[c][obj][thk] = thrust/weight * 100
                img_file = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.png'
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    #plot_superimposed_diagrams(form, form_base, save=img_file).show()
                    #plot_form(form, show_q=False, cracks=True).show()
            else:
                #plot_superimposed_diagrams(form, form_base).show()
                view = Viewer(form)
                view.show_solution()
                break

    view = Viewer(form, shape=vault)
    view.draw_thrust()
    view.draw_force()
    view.view_cracks()
    view.view_shape()
    view.view_reactions()
    view.show()

    import compas_tno
    location = compas_tno.get('form.json')
    form.to_json(location)

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])
