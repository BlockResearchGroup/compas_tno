from compas_plotters import Plotter
import matplotlib.pyplot as plt
from platform import node
from compas.geometry.intersections.intersections import intersection_segment_segment_xy
from compas.geometry import Line
import numpy as np
import math
# import cross_section_analysis as cra

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

import compas_tno

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import apply_bounds_tub_tlb

#  Parameters Geometry

type_formdiagram = 'arch'
type_structure = 'arch'

thk = 0.055
H = 0.68
L = 2.86
x0 = 0.0
B = 0.50  # m

x0_nobuttress = 0.60
xf_nobuttress = 2.26
h1_buttress = 0.51
h2_fill = 0.80
allow_thrust_buttress = False

add_load_mag = 2.0
span_percentage = 0.66

address_save = compas_tno.get('form.json')

discretisation_form = 20
discretisation_shape = 40

#  Parameters Reinforcement

fm = 14  # MPa
ffd = 390.08  # MPa
tf = 0.064
thk_mm = thk*1000  # mm
B_mm = B*1000  # mm

#  Parameters Construction  <- CHECK THIS IN THE PAPER

ro_bricks = 15.0
ro_butress = 12.0
ro_fill = 10.0

#  Parameters Optimisation

objective = 'max_section'  # try 'max'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 'tub']  # in the future add 'tlb' as variables
features = ['fixed']
starting_point = 'current'
tubmax = 0.5
tlbmax = 0.5

scalefactor = 0.5

data_shape = {
    'type': type_structure,
    'thk': thk,
    'H': H,
    'L': L,
    'x0': x0,
    'discretisation': discretisation_shape,
    'b': B,
    't': 0.0,
}

arch = Shape.from_library(data_shape)
arch.ro = ro_bricks
arch_swt = arch.compute_selfweight()

print('Arch Selfweight (kN):', arch_swt)

form = FormDiagram.create_arch(L=L, x0=x0, total_nodes=discretisation_form)

apply_envelope_from_shape(form, arch)
apply_selfweight_from_shape(form, arch)
apply_bounds_on_q(form, qmax=1e-6)
apply_bounds_tub_tlb(form, tubmax=tubmax, tlbmax=tlbmax)

# Compute and store vertex dx

for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    edges = form.vertex_edges(key)
    if len(edges) == 1:
        dx = abs(x - form.edge_midpoint(*edges[0])[0])
    else:
        dx = abs(form.edge_midpoint(*edges[0])[0] - form.edge_midpoint(*edges[1])[0])
    form.vertex_attribute(key, 'dx', dx)

# Add vertical load from buttress and filling

for key in form.vertices():
    x, _, _ = form.vertex_coordinates(key)
    zub = form.vertex_attribute(key, 'ub')
    dx = form.vertex_attribute(key, 'dx')
    print('\n', '-'*20)
    print('key:', key)
    print('zub:', zub)
    print('dx:', dx)
    print('x:', x)
    if (x < x0_nobuttress) or (x > xf_nobuttress):
        W_add_1 = dx * (h1_buttress - zub) * ro_butress * B
        W_add_2 = dx * (h2_fill - h1_buttress) * ro_fill * B
        if allow_thrust_buttress:
            form.vertex_attribute(key, 'ub', h1_buttress)
    else:
        W_add_1 = 0.0
        W_add_2 = dx * (h2_fill - zub) * ro_fill * B
    print('W_add_1:', W_add_1)
    print('W_add_2:', W_add_2)
    load_arch = form.vertex_attribute(key, 'pz')
    form.vertex_attribute(key, 'pz', load_arch - W_add_1 - W_add_2)
    print('load:', load_arch)
    print('load final:', form.vertex_attribute(key, 'pz'))

lumped_load = form.lumped_swt()
print('\nTotal Lumped Load is:', lumped_load)
print('Load from the filling is:', lumped_load - (- arch_swt))

# Find point @ 1/3 span and apply concentrated load.

distances_load = []
for key in form.vertices():
    x, _, _ = form.vertex_coordinates(key)
    dist = abs(x - span_percentage * L)
    distances_load.append(dist)

min_value = min(distances_load)
min_vertex = distances_load.index(min_value)

print('Vertex chosed to apply the load of {0} kN is vertex #{1}. This vertex is distant {2} to the point of application at {3} of the span'.format(add_load_mag, min_vertex, min_value, span_percentage))
original_load = form.vertex_attribute(min_vertex, 'pz')
form.vertex_attribute(min_vertex, 'pz', original_load - add_load_mag)
print('The old load in the vertex was: {0} | now the load is {1}'.format(original_load, original_load - add_load_mag))

# Run Optimisation

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['objective'] = objective
optimiser.settings['printout'] = True
optimiser.settings['find_inds'] = False
optimiser.settings['starting_point'] = starting_point

analysis = Analysis.from_elements(arch, form, optimiser)
analysis.set_up_optimiser()
analysis.run()

# Plot output

# form.to_json(address_save)
# print('Saved Solution to:', address_save)

view = Viewer(form, shape=arch)
view.draw_thrust()
view.draw_shape()
view.show()

view = Viewer(form)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
view.show()

# # ztest = get_shape_middle(arch, 2.5, 0)
# # print('ztest is:', ztest)

# # address = os.path.join(path, 'arch10.json')
# # arch = Shape.from_json(address)

# lines_normal = []
# points = []
# risultato_2 = []
# risultato_3 = []
# Aminis_2 = []
# Aminis_3 = []
# points2 = []
# node = []

# i = 0

# for vertex in arch.middle.vertices():
#     coord = arch.middle.vertex_coordinates(vertex)  # coordinates of the vertices of the thrust in the middle surface

#     # xi = arch.middle.vertex_attribute(u, 'x')
#     if coord[1] == 0:
#         norm_coord = arch.middle.vertex_normal(vertex)  # coordinates of the NORMAL to the vertices in the middle surface

#         p1 = [coord[0] - scalefactor*norm_coord[0], coord[2] - scalefactor*norm_coord[2]]  # starting point of the line that identifies the cross section
#         p2 = [coord[0] + scalefactor*norm_coord[0], coord[2] + scalefactor*norm_coord[2]]  # ending point of the line that identifies the cross section

#         # x1, z1 = coord[0], coord[2]
#         # x2, z2 = coord[0] - scalefactor*norm_coord[0], coord[2] - scalefactor*norm_coord[2]

#         l1 = Line(p1, p2)

#         for u, v in form.edges():
#             ucord = form.vertex_coordinates(u)
#             vcord = form.vertex_coordinates(v)
#             thrust_segment = [[ucord[0], ucord[2]], [vcord[0], vcord[2]]]

#             intersect = intersection_segment_segment_xy(l1, thrust_segment)

#             if intersect:

#                 alfirad = np.arctan((ucord[2]-vcord[2])/(ucord[0]-vcord[0]))  # inclination of the edge which intersect the cross-ection
#                 alfideg = abs(math.degrees(alfirad))
#                 betarad = np.arctan((p2[0]-p1[0])/(p2[1]-p1[1]))  # inclination of the cross-ection
#                 betadeg = abs(math.degrees(betarad))
#                 gammadeg = abs(alfideg - betadeg)  # deviation between the edge and the cross-section
#                 gammarad = abs(math.radians(gammadeg))

#                 force = -form.edge_attribute((u, v), 'f')*1000
#                 force_hor = force*np.cos(alfirad)  # normal force acting in the cross-section

#                 force_norm = force*np.cos(gammarad)  # normal force acting in the cross-section
#                 force_shear = force*np.sin(gammarad)  # shear force acting in the cross-section

#                 # tub = form.vertex_attribute(i, 'tub')*1000
#                 tlb = form.vertex_attribute(i, 'tlb')*1000

#                 # if tub > 10e-3:  #point out of the thickness -> need reinforcement
#                 #     moment_2 = force_hor * (tub + t/2)
#                 #     [Amini_2, wi] = cra.minimum_area(force_norm, moment_2, t, B, tub, fm, ffd, tf) #minimum area of reinforcement

#                 #     tb_norm = tub*np.cos(betarad)
#                 #     moment_3 = force_norm * (tb_norm + t/2)
#                 #     [Amini_3, wi] = cra.minimum_area(force_norm, moment_3, t, B, tb_norm, fm, ffd, tf) #minimum area of reinforcement

#                 if tlb > 10e-3:
#                     moment_2 = force_hor * (tlb + t/2)
#                     [Amini_2, wi] = cra.minimum_area(force_norm, moment_2, t, B, tlb, fm, ffd, tf)  # minimum area of reinforcement
#                     [Amini_2, wi] = [-Amini_2, wi]

#                     tb_norm = tlb*np.cos(betarad)
#                     moment_3 = force_norm * (tb_norm + t/2)
#                     [Amini_3, wi] = cra.minimum_area(force_norm, moment_3, t, B, tb_norm, fm, ffd, tf)  # minimum area of reinforcement
#                     [Amini_3, wi] = [-Amini_3, wi]

#                 else:  # point into the thickness -> NO need reinforcement
#                     [Amini_2, wi] = [0, 0]
#                     [Amini_3, wi] = [0, 0]

#                 risultato_2.append(
#                     {
#                         'normal node': vertex,
#                         'u':  u,
#                         'v': v,
#                         # 'tub' : round(tub,2),
#                         # 'tub_norm' : round(tub_norm,2),
#                         'force': str(round(force, 2)),
#                         'edge_angle': str(round(alfideg, 2)),
#                         'section_angle': str(round(betadeg, 2)),
#                         'gamma_angle': round(gammadeg, 2),
#                         'normal force': str(round(force_norm, 2)),
#                         'shear force': str(round(force_shear, 2)),
#                         'Amin_2': round(Amini_2, 2)
#                     }
#                 )

#                 risultato_3.append(
#                     {
#                         'normal node': vertex,
#                         'u':  u,
#                         'v': v,
#                         # 'tub' : round(tub,2),
#                         # 'tub_norm' : round(tub_norm,2),
#                         'force': str(round(force, 2)),
#                         'edge_angle': str(round(alfideg, 2)),
#                         'section_angle': str(round(betadeg, 2)),
#                         'gamma_angle': round(gammadeg, 2),
#                         'normal force': str(round(force_norm, 2)),
#                         'shear force': str(round(force_shear, 2)),
#                         'Amin_3': round(Amini_3, 2)
#                     }
#                 )

#                 node.append(i)
#                 Aminis_2.append(Amini_2)
#                 Aminis_3.append(Amini_3)

#                 points.append(
#                     {
#                         'pos': intersect,
#                         'radius': 0.02,
#                         # 'text': str(round(Amini,2))
#                     }
#                 )

#                 i = i+1

#                 break

#         lines_normal.append(
#             {
#                 'start': p1,
#                 'end': p2
#             }
#         )

# Area = plt.plot(node, Aminis_2, node, Aminis_3)
# plt.setp(Area, linewidth=2.0)
# plt.ylabel('reinforcement areas [mm^2]')
# plt.xlabel('cross-section', loc=None, fontsize=10)
# plt.xlim(0, 18)
# # plt.ylim(-5, 10)

# plt.show()

# # print(risultato)

# plotter = plot_form_xz(form, arch)
# # plotter.draw_lines(lines_normal)
# plotter.draw_points(points)
# # plotter.draw_points(points2)
# plotter.show()
