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

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis

from compas.geometry import normalize_vector
from numpy import array
import os

span = 10.0
k = 1.0
discretisation = 14
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 2 * discretisation
hc = None
he = None

save = True
solutions = {}

objective = ['Ecomp-linear']
objective = ['max', 'min']
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])
    lambd = 0.1

Xc = [5.0, 5.0, 0.0]
c = 0.1

for sign in [1]:  # get the minimum and maximum
    solutions[sign] = {}

    for obj in objective:  # set the objective
        solutions[sign][obj] = {}

        for thk in [0.5, 0.45, 0.40, 0.35, 0.33]:  # thickness of the problem

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

            solutions[sign][obj][thk] = {}

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

            for key in form.vertices_where({'is_fixed': True}):
                x, y, z = form.vertex_coordinates(key)
                dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), sign*(z - Xc[2])])
                vector_supports.append(dXbi)
                lines.append({
                    'start': [x, y, z],
                    'end': [x + dXbi[0], y + dXbi[1], z + dXbi[2]],
                    'width': 3
                })

            dXb = array(vector_supports)
            print(dXb)

            # from compas_plotters import MeshPlotter

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
            optimiser.settings['max_iter'] = 500
            optimiser.settings['gradient'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['printout'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['derivative_test'] = True
            optimiser.settings['starting_point'] = starting_point
            optimiser.settings['support_displacement'] = dXb

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(vault, form, optimiser)
            analysis.set_up_optimiser()
            analysis.run()

            form.overview_forces()
            if obj == 't':
                thk = form.attributes['thk']
            # if obj == 'Ecomp-linear':
            #     Ecomp = sign*optimiser.fopt

            weight = 0
            for key in form.vertices():
                weight += form.vertex_attribute(key, 'pz')

            thrust = form.thrust()

            T_over_W = abs(thrust/weight)
            # E_over_W = abs(Ecomp/weight)

            print('Ratio Thrust/Weight:', T_over_W)
            # print('Ratio Ecomp/Weight:', E_over_W)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'compl_energy', type_structure, type_formdiagram)
            if 'fixed' in optimiser.settings['features']:
                folder = os.path.join(folder, 'fixed')
            else:
                folder = os.path.join(folder, 'mov_c_' + str(c))
            if he:
                folder = os.path.join(folder, 'hc_' + str(hc) + '_he_' + str(he))
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
                solutions[sign][obj][thk]['T/W'] = T_over_W
                # solutions[sign][obj][thk]['E/W'] = E_over_W
                if save:
                    form.to_json(address)
                    print('Saved to: ', address)
                    plot_superimposed_diagrams(form, form_base, save=img_file)#.show()
                    plot_form(form, show_q=False, cracks=True)#.show()
            else:
                # plot_superimposed_diagrams(form, form_base).show()
                # view = Viewer(form)
                # view.show_solution()
                break

view = Viewer(form)
view.show_solution()

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])
