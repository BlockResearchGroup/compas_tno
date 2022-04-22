from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.plotters import TNOPlotter
from compas.geometry import Vector
from compas.geometry import Point

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from compas.geometry import normalize_vector
from numpy import array
import math
import os

span = 10.0
k = 1.0
discretisation = 50
type_formdiagram = 'arch'  # write the type of form diagram you want and is in the file shape
type_structure = 'arch'
thk = 1.0

save = True
solutions = {}

objective = ['Ecomp-linear']
Ecomp_method = 'complete'
solver = 'SLSQP'
constraints = ['funicular', 'envelope', 'reac_bounds']  # , 'envelopexy'
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None
starting_point = 'current'

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambdh']:
    variables.append(objective[0])
    lambd = 0.1

Xc = [0.5, 0.0, 0.0]

c = 0

for i_angle in range(36):  # set the distance that the nodes can move
    theta = math.radians(i_angle*10.0)
    solutions[c] = {}

    for obj in objective:  # set the objective
        solutions[c][obj] = {}

        for thk in [thk]:  # thickness of the problem

            # Create form diagram

            data_diagram = {
                'type': type_formdiagram,
                'H': span/2,
                'L': span,
                'x0': 0,
                'total_nodes': discretisation,
            }

            form = FormDiagram.from_library(data_diagram)

            # Create shape

            data_shape = {
                'type': type_structure,
                'thk': thk,
                'H': span/2,
                'L': span,
                'x0': 0,
                'discretisation': 100,
                'b': 1.0,
                't': 0.5,
            }

            shape = Shape.from_library(data_shape)
            shape.ro = 20.0

            # ------------------------------------------------------------
            # -----------------------  INITIALISE   ----------------------
            # ------------------------------------------------------------

            # Apply Selfweight and Envelope

            apply_bounds_on_q(form, qmax=0.0)

            form_base = form.copy()

            # Consider the displacement vector

            if 'Ecomp' in obj.split('-'):

                lines = []
                vector_supports = []
                vectors_plot = []
                base_plot = []

                sign = None

                for key in form.vertices_where({'is_fixed': True}):
                    x, y, z = form.vertex_coordinates(key)

                    if x < Xc[0]:
                        dXbi = [0, 0, 0]  # left support of the arch
                    else:
                        dXbi = [math.cos(theta), 0, math.sin(theta)]

                    vector_supports.append(dXbi)
                    # vectors_plot.append(Vector(dXbi[0], dXbi[1], dXbi[2]))
                    vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
                    base_plot.append(Point(x, y, z - 0.2))
                    lines.append({
                        'start': [x, y, z],
                        'end': [x + dXbi[0], y + dXbi[1], z + dXbi[2]],
                        'width': 3
                    })

                dXb = array(vector_supports)
                print(dXb)

                key_index = form.key_index()
                # plotter = TNOPlotter(form, shape=shape)
                # plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
                # plotter.draw_shape()
                # plotter.draw_vectors(vectors=vectors_plot, bases=base_plot)
                # # plotter.draw_vertices(keys=form.fixed(), facecolor={key: '000000' for key in form.vertices_where({'is_fixed': True})})
                # # plotter.draw_arrows(lines)
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
            optimiser.settings['axis_sym'] = axis_sym
            optimiser.settings['max_iter'] = 500
            optimiser.settings['gradient'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['printout'] = True
            optimiser.settings['jacobian'] = True
            optimiser.settings['derivative_test'] = False
            optimiser.settings['starting_point'] = starting_point
            optimiser.settings['support_displacement'] = dXb
            optimiser.settings['Ecomp_method'] = Ecomp_method

            # --------------------- 5. Set up and run analysis ---------------------

            analysis = Analysis.from_elements(shape, form, optimiser)
            analysis.apply_selfweight()
            analysis.apply_envelope()
            analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            analysis.run()

            form.overview_forces()
            if obj == 't':
                thk = form.attributes['thk']

            weight = 0
            for key in form.vertices():
                weight += form.vertex_attribute(key, 'pz')

            CEnergy = 0.0

            for i, key in enumerate(form.vertices_where({'is_fixed': True})):
                rx = form.vertex_attribute(key, '_rx')
                ry = form.vertex_attribute(key, '_ry')
                rz = form.vertex_attribute(key, '_rz')
                CEnergy += - dXb[i][0] * rx - dXb[i][1] * ry - dXb[i][2] * rz

            thrust = form.thrust()
            print('Thrust/2:', thrust/2)
            print('Thrust/Weight/2:', thrust/weight/2)

            print('CEnergy:', CEnergy)
            print('Ratio CEnergy/Weight:', CEnergy/weight)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'compl_energy', 'assym', type_structure, type_formdiagram)
            if 'ind' in optimiser.settings['variables']:
                folder = os.path.join(folder, 'fixed')
            else:
                if 'fixed' not in features:
                    folder = os.path.join(folder, 'mov_c_' + str(c))
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)
            address = save_form + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'

            print('Optimiser exitflag:', optimiser.exitflag)

            p0 = [-1.1, -1.1]
            p1 = [11.1, -1.1]
            p2 = [11.1, 6]
            p3 = [-1.1, 6]
            lines_around = [[p0, p1], [p1, p2], [p2, p3], [p3, p0]]

            folder = os.path.join(folder, 'img')
            pic = folder + '/fig-' + str(i_angle) + '.png'
            print('Picture saved to:', pic)

            plotter = TNOPlotter(form, shape=shape, figsize=(12, 8))
            plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
            plotter.settings['show.reactions.extended'] = True
            plotter.settings['show.reactions.asarrows'] = False
            plotter.draw_form_xz()
            plotter.draw_vectors(vectors=vectors_plot, bases=base_plot)
            plotter.draw_shape_xz()
            plotter.draw_cracks()
            plotter.draw_reactions()
            plotter.draw_lines(lines=lines_around)
            plotter.save(pic)
            # plotter.show()

# plotter = TNOPlotter(form, shape=shape)
# plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
# plotter.settings['show.reactions.extended'] = True
# plotter.settings['show.reactions.asarrows'] = False
# plotter.draw_form_xz()
# plotter.draw_vectors(vectors=vectors_plot, bases=base_plot)
# plotter.draw_shape_xz(stereotomy=True, blocks=20)
# plotter.draw_cracks()
# plotter.draw_reactions()
# plotter.show()


plotter = TNOPlotter(form, shape=shape)
plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
plotter.settings['show.reactions.extended'] = True
plotter.settings['show.reactions.asarrows'] = False
plotter.draw_form_xz()
plotter.draw_vectors(vectors=vectors_plot, bases=base_plot)
plotter.draw_shape_xz()
plotter.draw_cracks()
plotter.draw_reactions()
plotter.draw_lines(lines=lines_around)
plotter.show()

print(solutions)
print('\n')

for key in solutions:
    for key2 in solutions[key]:
        print(key, key2, solutions[key][key2])
