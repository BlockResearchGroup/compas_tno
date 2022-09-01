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
thk = 1.0  # t/R = 0.20
thk = 0.75  # t/R = 0.15

save = False
solutions = {}

objective = ['Ecomp-linear']
Ecomp_method = 'complete'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']  # , 'envelopexy'
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None
starting_point = 'loadpath'

Xc = [0.5, 0.0, 0.0]

c = 0

# for i_angle in [0, 9, 18, 27]: # range(36):  # set the distance that the nodes can move
for i_angle in [11, 25]: # range(36):  # set the distance that the nodes can move
# for i_angle in range(36):  # range(36):  # set the distance that the nodes can move
    theta = math.radians(i_angle*10.0)
    solutions[i_angle] = {}

    for obj in objective:  # set the objective
        # solutions[i_angle][obj] = {}

        for thk in [thk]:  # thickness of the problem

            form = FormDiagram.create_arch(H=span/2, L=span, x0=0.0, discretisation=discretisation)

            shape = Shape.create_arch(H=span/2, L=span, thk=thk, x0=0.0, t=0.5)

            if 'Ecomp' in obj.split('-'):

                vector_supports = []
                vectors_plot = []
                base_plot = []
                sup_key = None

                sign = None

                for key in form.vertices_where({'is_fixed': True}):
                    x, y, z = form.vertex_coordinates(key)

                    if x < Xc[0]:
                        dXbi = [0, 0, 0]  # left support of the arch
                        sup_key = key
                    else:
                        dXbi = [math.cos(theta), 0, math.sin(theta)]
                        sup_key = key

                    vector_supports.append(dXbi)
                    vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
                    base_plot.append(Point(x, y, z - 0.2))

                dXb = array(vector_supports)
                print(dXb)

                key_index = form.key_index()

            analysis = Analysis.create_compl_energy_analysis(form,
                                                             shape,
                                                             printout=True,
                                                             support_displacement=dXb,
                                                             max_iter=10000,
                                                             Emethod=Ecomp_method,
                                                             solver=solver,
                                                             starting_point=starting_point)

            analysis.optimiser.set_constraints(constraints)
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

            print('Weight:', weight)
            left_support = []
            right_support = []

            CEnergy = 0.0

            left = True
            for i, key in enumerate(form.vertices_where({'is_fixed': True})):
                rx = form.vertex_attribute(key, '_rx')
                ry = form.vertex_attribute(key, '_ry')
                rz = form.vertex_attribute(key, '_rz')
                x, y, z = form.vertex_coordinates(key)
                if left:
                    left_support = [x, rx, ry, rz]
                    left = False
                else:
                    right_support = [x, rx, ry, rz]
                    T = rx
                    V = rz
                CEnergy += - dXb[i][0] * rx - dXb[i][1] * ry - dXb[i][2] * rz

            print('left', left_support)
            print('right', right_support)

            solutions[i_angle]['T/W'] = rx
            solutions[i_angle]['V/W'] = rz
            solutions[i_angle]['weight'] = abs(weight)
            print(vector_supports)
            solutions[i_angle]['ux'] = vector_supports[1][0]
            solutions[i_angle]['uz'] = vector_supports[1][2]

            # solutions[i_angle]['CEnergy'] = CEnergy
            # solutions[i_angle]['lp'] = form.loadpath()

            # thrust = form.thrust()
            # print('Thrust/2:', thrust/2)
            # print('Thrust/Weight/2:', thrust/weight/2)

            # print('CEnergy:', CEnergy)
            # print('Ratio CEnergy/Weight:', CEnergy/weight)

            # folder = os.path.join('/Users/mricardo/compas_dev/me', 'compl_energy', 'assym', type_structure, type_formdiagram)
            # if 'ind' in optimiser.settings['variables']:
            #     folder = os.path.join(folder, 'fixed')
            # os.makedirs(folder, exist_ok=True)
            # title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            # save_form = os.path.join(folder, title)
            # address = save_form + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'

            # print('Optimiser exitflag:', analysis.optimiser.exitflag)

            margin = 2.0
            p0 = [-margin, -margin]
            p1 = [span + margin, -margin]
            p2 = [span + margin, span/2 + margin]
            p3 = [-margin, span/2 + margin]
            lines_around = [[p0, p1], [p1, p2], [p2, p3], [p3, p0]]

            for u, v in form.edges():
                uz = form.vertex_attribute(u, 'z')
                vz = form.vertex_attribute(v, 'z')

            print('Plot for i_angle', i_angle)

            if save:
                path = '/Users/mricardo/compas_dev/me/compl_energy/arch/solution-i={}.json'.format(i_angle)
                form.to_json(path)
                print('Saved to:', path)

            # folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/arch/imgs')
            # pic = folder + '/fig-' + str(i_angle) + '.png'
            # print('Picture saved to:', pic)

            # plotter = TNOPlotter(form, shape=shape, figsize=(10, 8))
            # # plotter = TNOPlotter(form, shape=shape, figsize=(8, 4))
            # plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
            # plotter.settings['show.reactions.extended'] = True
            # plotter.settings['show.reactions.asarrows'] = False
            # plotter.settings['size.edge.base_thickness'] = 2.0
            # plotter.settings['scale.reactions'] = 0.02
            # plotter.draw_form_xz(scale_width=False)
            # for i in range(len(vectors_plot)):
            #     vector = vectors_plot[i]
            #     base = base_plot[i]
            #     plotter.draw_vector(vector=vector, base=base)
            # plotter.draw_shape_xz()
            # plotter.draw_cracks()
            # plotter.draw_reactions(scale_width=False)
            # plotter.draw_lines(lines=lines_around)
            # # plotter.save(pic)
            # plotter.show()


for i in solutions:
    print(i, solutions[i]['weight'], solutions[i]['T/W'], solutions[i]['V/W'], solutions[i]['ux'], solutions[i]['uz'])

