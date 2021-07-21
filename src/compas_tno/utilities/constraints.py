from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from compas.utilities import geometric_key
from compas.utilities import reverse_geometric_key

from compas.geometry import closest_point_on_line

from compas_tno.algorithms.equilibrium import z_from_form

from compas_plotters import MeshPlotter
from copy import deepcopy

import math

from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_form_xz


__all__ = [
    'check_constraints',
    'distance_target',
    'replicate_contraints',
    'interp_surf',
    'null_edges',
    'rectangular_smoothing_constraints',
]


def check_constraints(form, show=False, lb_show=False, ub_show=False, tol=1e-6):

    try:
        t = form.attributes['offset']
    except:
        t = 0.0
    outside = {}
    lbs = {}
    ubs = {}
    penalty = 0

    for key, vertex in form.vertex.items():
        z = form.vertex_coordinates(key)[2] + t
        if vertex.get('lb', None):
            lb = vertex['lb']
            lbs[key] = lb
            if z < lb - tol:
                outside[key] = lb - z
                penalty += (abs(outside[key])+4)**(4)
        if vertex.get('ub', None):
            ub = vertex['ub']
            ubs[key] = ub
            if z > ub + tol:
                outside[key] = z - ub
                penalty += (abs(outside[key])+4)**(4)

    print('The penalty in the constraints is {0:.3f}'.format(penalty))

    if penalty > 0.1:
        if show:
            plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
            plotter.draw_vertices(text=outside)
            plotter.draw_edges()
            plotter.show()

        if lb_show:
            plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
            plotter.draw_vertices(text=lbs)
            plotter.draw_edges()
            plotter.show()

        if ub_show:
            plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
            plotter.draw_vertices(text=ubs)
            plotter.draw_edges()
            plotter.show()

    return penalty


def distance_target(form, method='least-squares'):

    dist = 0

    if method == 'least-squares':

        for key in form.vertices():
            _, _, z = form.vertex_coordinates(key)
            targ = form.vertex_attribute(key, 'target')
            dist += (z - targ) ** 2

    if method == 'volume':

        for key in form.vertices():
            _, _, z = form.vertex_coordinates(key)
            targ = form.vertex_attribute(key, 'target')
            weight = form.vertex_attribute(key, 'weight')
            dist += abs(z - targ)*weight

    return dist


def replicate_contraints(file, file_constraint):

    form = FormDiagram.from_json(file)
    form_ = FormDiagram.from_json(file_constraint)
    gkey_planar = {}

    for key_real in form.vertices():
        coord = form.vertex_coordinates(key_real)
        gkey_proj = geometric_key([coord[0], coord[1], 0.0])
        gkey_planar[gkey_proj] = key_real

    for key in form_.vertices():
        target = form_.vertex[key].get('target', 0.0)
        lb = form_.vertex[key].get('lb', 0.0)
        ub = form_.vertex[key].get('ub', 0.0)
        if target < 10**(-4):
            target = 0.00
        gkey = geometric_key([form_.vertex_coordinates(key)[0], form_.vertex_coordinates(key)[1], 0.0])
        form.vertex_attribute(gkey_planar[gkey], 'target', target)
        form.vertex_attribute(gkey_planar[gkey], 'lb', lb)
        form.vertex_attribute(gkey_planar[gkey], 'ub', ub)

    return form


def interp_surf(form):

    x = []
    y = []
    s = []

    for key, vertex in form.vertex.items():
        if vertex.get('_is_external') == False:
            x.append(vertex.get('x'))
            y.append(vertex.get('y'))
            s.append(vertex.get('target'))

    from scipy import interpolate

    surf = interpolate.interp2d(x, y, s, kind='linear')

    return surf


def null_edges(form, plot=False):

    null_edges = []
    all_edges = []

    for u, v in form.edges():
        if form.edge_attribute((u, v), '_is_external') == False and form.edge_attribute((u, v), '_is_edge') == True:
            activ = 0
            coord_u = form.vertex_coordinates(u)
            coord_v = form.vertex_coordinates(v)
            ux = round(coord_u[0], 3)
            uy = round(coord_u[1], 3)
            vx = round(coord_v[0], 3)
            vy = round(coord_v[1], 3)
            mid_x, mid_y, _ = form.edge_midpoint(u, v)
            if uy == vy and ((uy is not 10.0 and vy is not 10.0) or (uy is not 0.0 and vy is not 0.0)):
                if (mid_y > mid_x and mid_y < 10 - mid_x) or (mid_y < mid_x and mid_y > 10 - mid_x):
                    if uy == 5.0 and vy == 5.0 and ux > 0.01 and vx > 0.01 and ux < 9.99 and vx < 9.99:  # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u, v))
                        activ += 1
            if ux == vx and ((ux is not 10.0 and vx is not 10.0) or (ux is not 0.0 and vx is not 0.0)):
                if (mid_y > mid_x and mid_y > 10 - mid_x) or (mid_y < mid_x and mid_y < 10 - mid_x):
                    if ux == 5.0 and vx == 5.0 and uy > 0.01 and vy > 0.01 and uy < 9.99 and vy < 9.99:  # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u, v))
                        activ += 1
            if activ == 0:
                all_edges.append((u, v))

    if plot:
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges(all_edges)
        plotter.draw_edges(null_edges, color='ff0000')
        plotter.show()

    return null_edges


def rectangular_smoothing_constraints(form, xy_span=[[0, 10], [0, 10]]):

    [[x0, x1], [y0, y1]] = xy_span
    cons = {key: None for key in form.vertices()}
    line_top = [[x0, y1, 0], [x1, y1, 0]]
    line_bottom = [[x0, y0, 0], [x1, y0, 0]]
    line_left = [[x0, y0, 0], [x0, y1, 0]]
    line_right = [[x1, y0, 0], [x1, y1, 0]]
    for key in form.vertices_on_boundary():
        x, y, z = form.vertex_coordinates(key)
        if x == x0:
            cons[key] = line_left
        elif x == x1:
            cons[key] = line_right
        elif y == y0:
            cons[key] = line_bottom
        elif y == y1:
            cons[key] = line_top
    for key in form.vertices_where({'is_fixed': True}):
        cons[key] = form.vertex_coordinates(key)
    return cons


def create_cracks(form, dx=[[0.50, 0.55]], dy=[[-0.1, 0.1]], type=['top'], view=False):
    """ Create cracks on a form diagram to the nodes desired.

    Parameters
    ----------
    form : obj
        ForceDiagram to constraint.
    dx : list
        List of the range on the x-coordinates of the nodes to be constrained.
    dy : list
        List of the range on the y-coordinates of the nodes to be constrained.
    type : list
        List with the type of constraint to applye to the vertices (top or bottom).

    Returns
    -------
    form : obj
        ForceDiagram with the constraints in attribute 'cracks'.

    """

    k_i = form.key_index()
    cracks_ub = []
    cracks_lb = []

    count = 0

    for key in form.vertices():
        x, y, _ = form.vertex_coordinates(key)
        for i in range(len(dx)):
            if dx[i][0] <= x <= dx[i][1] and dy[i][0] <= y <= dy[i][1]:
                if type[i] == 'top':
                    cracks_ub.append(k_i[key])
                if type[i] == 'bottom':
                    cracks_lb.append(k_i[key])

    form.attributes['cracks'] = (cracks_lb, cracks_ub)

    if view:
        plot_form_xz(form, radius=0.02, cracks=True).show()

    return form


def circular_joints(form, x0=None, xf=None, blocks=18, thk=0.5, t=0.0, tol=1e-3):

    k_i = form.key_index()

    if x0 == None or xf == None:
        x = []
        for key in form.vertices():
            x.append(form.vertex_coordinates(key)[0])
        x0 = min(x)
        xf = max(x)
    y = 0.0

    xc = (xf+x0)/2
    r = xf - xc
    ri = r - thk/2
    re = r + thk/2
    form.attributes['Re'] = re
    form.attributes['Ri'] = ri
    print('SpanMid: {0:.2} m / SpanInt: {1:.2} m / SpanExt: {2:.2} m / Thickness: {3:.4} m / Ratio t/Ri: {4:.4} m / Ratio t/R: {5:.4} m / Number of Blocks: {6}'.format(2 *
          r, 2*ri, 2*re, thk, (thk/ri), (thk/r), blocks))

    njoints = blocks+1
    joints = {}
    for j in range(njoints):
        theta = j/blocks*math.pi
        xi = xc + ri * math.cos(theta)  # takeout
        zi = ri * math.sin(theta) + tol
        xe = xc + re * math.cos(theta)
        ze = re * math.sin(theta) + tol  # take out
        xmax = max(xi, xe)
        xmin = min(xi, xe)
        possible_edges = []
        for (u, v) in form.edges():
            xu, xv = form.vertex_coordinates(u)[0], form.vertex_coordinates(v)[0]
            if max(xu, xv) >= xmin and min(xu, xv) <= xmax:
                possible_edges.append(tuple(sorted([k_i[u], k_i[v]])))
                if form.vertex_attribute(u, 'is_fixed') == True:
                    possible_edges.append(tuple(sorted([-k_i[u], k_i[u]])))
                if form.vertex_attribute(v, 'is_fixed') == True:
                    possible_edges.append(tuple(sorted([-k_i[v], k_i[v]])))
        joints[j] = [[xi, y, zi], [xe, y, ze], set(possible_edges)]
        print(joints[j])
    form.attributes['joints'] = joints

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        zt = math.sqrt(r**2 - (x-xc)**2)
        ze = math.sqrt(re**2 - (x-xc)**2) - t
        form.vertex_attribute(key, 'target', value=zt)
        zi2 = ri**2 - (x-xc)**2
        if zi2 < 0:
            zi = 0 - t
        else:
            zi = math.sqrt(zi2) - t
        if form.vertex_attribute(key, 'is_fixed') == True:
            form.vertex_attribute(key, 'lb', value=None)
            form.vertex_attribute(key, 'ub', value=None)
        else:
            form.vertex_attribute(key, 'lb', value=zi)
            form.vertex_attribute(key, 'ub', value=ze)
        # form.vertex_attribute(key,'z',value=ze)
        if form.vertex_attribute(key, 'is_fixed') == True:
            form.vertex_attribute(key, 'b', value=[thk/2, 0.0])
        if x == x0:
            form.attributes['tmax'] = ze

    return form


def rollers_on_openings(form, xy_span=[[0.0, 10.0], [0.0, 10.0]], max_f=5.0, constraint_directions='all'):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    bndr = form.vertices_on_boundary()

    for key in bndr:
        if form.vertex_attribute(key, 'is_fixed') == False:
            x, y, _ = form.vertex_coordinates(key)
            if x == x1 and (constraint_directions in ['all', 'x']):
                form.vertex_attribute(key, 'rol_x', True)
                form.vertex_attribute(key, 'max_rx', max_f)
            if x == x0 and (constraint_directions in ['all', 'x']):
                form.vertex_attribute(key, 'rol_x', True)
                form.vertex_attribute(key, 'max_rx', max_f)
            if y == y1 and (constraint_directions in ['all', 'y']):
                form.vertex_attribute(key, 'rol_y', True)
                form.vertex_attribute(key, 'max_ry', max_f)
            if y == y0 and (constraint_directions in ['all', 'y']):
                form.vertex_attribute(key, 'rol_y', True)
                form.vertex_attribute(key, 'max_ry', max_f)

    return form
