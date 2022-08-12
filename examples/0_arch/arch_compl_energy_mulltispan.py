from re import L

from cvxpy import Solution
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.utilities.envelopes import apply_bounds_reactions
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q
from compas.geometry import mirror_points_line

from compas.geometry import norm_vector
from compas.geometry import Line

from compas_tno.plotters import TNOPlotter
from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Line
from compas_view2.objects import Arrow

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from compas.geometry import normalize_vector
from compas.datastructures import Mesh
from compas.colors import Color
from numpy import array
import math
import os


def draw_multispan_lines(plotter, modules=2, H=1.00, L=2.0, x0=0.0, thk=0.20, total_nodes=50, stereotomy=False, close_bottom=True):
    """Helper to draw the lines of intrados and extrados of an arch for given parameters

    Parameters
    ----------
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0
    thk : float, optional
        Thickness of the arch, by default 0.20
    total_nodes : int, optional
        Density of the shape, equals to the number of blocks, by default 50
    stereotomy : bool, optional
        Whether or not interfaces of the stereotomy should be drawn, by default False
    close_bottom : bool, optional
        Whether or not the last interfaces of the arch should be drawn, by default True

    Returns
    -------
    None
        Lines are added to the plotter
    """

    lines_intrados = []
    lines_extrados = []
    line_sym = [[L, 0, 0], [L, H, 0]]

    radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2
    spr = math.atan2((L/2), (radius - H))
    tot_angle = 2*spr
    an_intra = tot_angle / (total_nodes)  # semicircular only

    theta = math.acos(radius/re)
    angle_init = (math.pi - tot_angle)/2
    angle_end = math.pi - theta
    tot_angle = (angle_end - angle_init)
    an_extra = tot_angle / (total_nodes)
    zc = radius - H
    xc = L/2 + x0
    i = 0

    for i in range(total_nodes):
        angle_ii = angle_init + i * an_intra
        angle_fi = angle_init + (i + 1) * an_intra
        angle_ie = angle_init + i * an_extra
        angle_fe = angle_init + (i + 1) * an_extra

        xii = xc - ri * math.cos(angle_ii)
        xif = xc - ri * math.cos(angle_fi)
        zii = ri * math.sin(angle_ii) - zc
        zif = ri * math.sin(angle_fi) - zc

        xei = xc - re * math.cos(angle_ie)
        xef = xc - re * math.cos(angle_fe)
        zei = re * math.sin(angle_ie) - zc
        zef = re * math.sin(angle_fe) - zc

        line_i = Line([xii, zii, 0.0], [xif, zif, 0.0])
        lines_intrados.append(line_i)

        line_e = Line([xei, zei, 0.0], [xef, zef, 0.0])
        lines_extrados.append(line_e)

    interfaces = []

    if close_bottom:
        interfaces = [0, total_nodes]
    if stereotomy:
        interfaces = list(range(total_nodes + 1))

    for i in interfaces:
        angle = angle_init + i * an_intra

        xi = xc - ri * math.cos(angle)
        xf = xc - re * math.cos(angle)
        zi = ri * math.sin(angle) - zc
        zf = re * math.sin(angle) - zc

        line = Line([xi, zi, 0.0], [xf, zf, 0.0])
        lines_intrados.append(line)

    # update this when collections are available

    for le in lines_extrados:
        plotter.app.add(le, draw_as_segment=True, color=Color.black())
        lesym = Line(*mirror_points_line(le, line=line_sym))
        plotter.app.add(lesym, draw_as_segment=True, color=Color.black())
    for li in lines_intrados:
        plotter.app.add(li, draw_as_segment=True, color=Color.black())
        lisym = Line(*mirror_points_line(li, line=line_sym))
        plotter.app.add(lisym, draw_as_segment=True, color=Color.black())

    return


span = 10.0
k = 1.0
discretisation = 20
type_formdiagram = 'arch'
type_structure = 'arch'
thk = 1.0

x0_1 = 0.0
x0_2 = 10.0

form1 = FormDiagram.create_arch(H=span/2, L=span, x0=x0_1, discretisation=discretisation)
arch1 = Shape.create_arch(H=span/2, L=span, thk=thk, x0=x0_1, t=0.5)

form2 = FormDiagram.create_arch(H=span/2, L=span, x0=x0_2, discretisation=discretisation)
arch2 = Shape.create_arch(H=span/2, L=span, thk=thk, x0=x0_2, t=0.5)

apply_selfweight_from_shape(form1, arch1)
apply_envelope_from_shape(form1, arch1)
apply_bounds_reactions(form1, arch1)

apply_selfweight_from_shape(form2, arch2)
apply_envelope_from_shape(form2, arch2)
apply_bounds_reactions(form2, arch2)

lines = []

for edge in form1.edges():
    lines.append(form1.edge_coordinates(*edge))
for edge in form2.edges():
    lines.append(form2.edge_coordinates(*edge))

mesh = Mesh.from_lines(lines)
form : FormDiagram = FormDiagram.from_mesh(mesh)

ubpts = []
lbpts = []
base = None
for key in form.vertices():
    if key < discretisation:
        keybase = key
        base = form1
    else:
        keybase = key - discretisation + 1
        base = form2

    ub = base.vertex_attribute(keybase, 'ub')
    b = base.vertex_attribute(keybase, 'b')
    lb = base.vertex_attribute(keybase, 'lb')
    target = base.vertex_attribute(keybase, 'target')
    pz = base.vertex_attribute(keybase, 'pz')
    fixed = base.vertex_attribute(keybase, 'is_fixed')

    form.vertex_attribute(key, 'ub', ub)
    form.vertex_attribute(key, 'lb', lb)
    form.vertex_attribute(key, 'target', target)
    form.vertex_attribute(key, 'pz', pz)
    form.vertex_attribute(key, 'b', b)
    form.vertex_attribute(key, 'is_fixed', fixed)

    ubpts.append(Point(form.vertex_attribute(key, 'x'), form.vertex_attribute(key, 'ub')))
    lbpts.append(Point(form.vertex_attribute(key, 'x'), form.vertex_attribute(key, 'lb')))

bridge = Shape.from_formdiagram_and_attributes(form)

# Plots

# plotter = TNOPlotter(form, shape=bridge, figsize=(12, 8))
# plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
# plotter.settings['show.reactions.extended'] = True
# plotter.settings['show.reactions.asarrows'] = False
# plotter.draw_form_xz()
# plotter.draw_supports()
# plotter.draw_shape_xz()
# for point in ubpts:
#     plotter.add(point, facecolor=Color.blue())
# for point in lbpts:
#     plotter.add(point, facecolor=Color.green())
# plotter.show()

solutions = {}
sup_key = None

for i_angle in range(36):  # set the distance that the nodes can move
    theta = math.radians(i_angle*10.0)
    solutions[i_angle] = {}

    xc = span
    vector_supports = []
    lines = []
    plots_vectors = []
    for key in form.vertices_where({'is_fixed': True}):

        x, y, z = form.vertex_coordinates(key)
        dXbi = [0, 0, 0]
        # if abs(x - xc) < 1e-3:
        if x - xc > 1:
            # dXbi = [0, 0, -1]
            dXbi = [math.cos(theta), 0, math.sin(theta)]
            sup_key = key
        if norm_vector(dXbi) > 1e-3:
            start = [x, y, z]
            end = [x + dXbi[0], y + dXbi[1], z + dXbi[2]]
            lines.append(Line(start, end))
            plots_vectors.append([Point(*start), Vector(dXbi[0], dXbi[2])])
        vector_supports.append(dXbi)

    dXb = array(vector_supports)
    print(dXb)

    # plotter : TNOPlotter = TNOPlotter(form, figsize=(12, 8))
    # # plotter.draw_form_xz()
    # draw_multispan_lines(plotter, H=span/2, thk=thk, L=span, x0=x0_1)
    # plotter.show()

    analysis  = Analysis.create_compl_energy_analysis(form, bridge, printout=True, support_displacement=dXb, max_iter=2000, solver='SLSQP')

    analysis.optimiser.set_constraints(['envelope', 'funicular', 'reac_bounds'])
    analysis.set_up_optimiser()
    analysis.run()

    weight = 0
    for key in form.vertices():
        weight += form.vertex_attribute(key, 'pz')

    CEnergy = 0.0

    for i, key in enumerate(form.vertices_where({'is_fixed': True})):
        rx = form.vertex_attribute(key, '_rx')
        ry = form.vertex_attribute(key, '_ry')
        rz = form.vertex_attribute(key, '_rz')
        CEnergy += - dXb[i][0] * rx - dXb[i][1] * ry - dXb[i][2] * rz

    T = form.vertex_attribute(sup_key, '_rx')
    V = form.vertex_attribute(sup_key, '_rz')

    solutions[i_angle]['T/W'] = abs(T/weight)
    solutions[i_angle]['V/W'] = abs(V/weight)
    solutions[i_angle]['weight'] = abs(weight)
    solutions[i_angle]['CEnergy'] = CEnergy
    solutions[i_angle]['lp'] = form.loadpath()

    p0 = [-1.1, -1.1]
    p1 = [21.1, -1.1]
    p2 = [21.1, 6]
    p3 = [-1.1, 6]
    lines_around = [[p0, p1], [p1, p2], [p2, p3], [p3, p0]]

    # print(solutions)
    print('Plot for i_angle', i_angle)

    folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/multispan/imgs')
    pic = folder + '/fig-' + str(i_angle) + '.png'
    print('Picture saved to:', pic)

    plotter : TNOPlotter = TNOPlotter(form, shape=bridge, figsize=(12, 8))
    plotter.settings['color.edges.shape'] = (0.0, 0.0, 0.0)
    plotter.settings['show.reactions.extended'] = True
    plotter.settings['show.reactions.asarrows'] = False
    plotter.draw_form_xz()
    # plotter.draw_supports()
    draw_multispan_lines(plotter, H=span/2, thk=thk, L=span, x0=x0_1)
    # plotter.draw_shape_xz()
    plotter.draw_reactions()
    for base, vector in plots_vectors:
        plotter.draw_vector(vector, base)
    plotter.draw_cracks()
    plotter.draw_lines(lines=lines_around)
    # plotter.save(pic)
    plotter.show()


for i in solutions:
    print(i, solutions[i]['weight'], solutions[i]['CEnergy'], solutions[i]['T/W'], solutions[i]['V/W'], solutions[i]['lp'])
