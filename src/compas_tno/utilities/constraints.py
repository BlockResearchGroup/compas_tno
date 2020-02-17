from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas.utilities import geometric_key
from compas.utilities import reverse_geometric_key

from compas.geometry.transformations.transformations import mirror_point_line
from compas.geometry.transformations.transformations import rotate_points
from compas.geometry import closest_point_on_line
from compas.geometry import midpoint_line_xy
from compas.geometry import matrix_from_axis_and_angle

from compas.geometry.distance import distance_point_point_xy
from numpy import argmin
from compas_tno.algorithms.equilibrium import z_from_form

from compas_plotters import MeshPlotter
from copy import deepcopy

import math

from compas_tno.plotters.plotters import plot_form
from compas_tno.plotters.plotters import plot_form_xz

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'check_constraints',
    'distance_target',
    'replicate_contraints',
    'interp_surf',
    'null_edges',
    'set_height_constraint',
    'set_cross_vault_heights',
    'set_pavillion_vault_heights',
    'set_oct_vault_heights',
    'set_dome_heights',
    'circular_heights',
    'circular_joints',
]


def check_constraints(form, show=False, lb_show=False, ub_show=False, tol = 1e-6):

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
            targ = form.get_vertex_attribute(key, 'target')
            dist += (z - targ) ** 2

    if method == 'volume':

        for key in form.vertices():
            _, _, z = form.vertex_coordinates(key)
            targ = form.get_vertex_attribute(key, 'target')
            weight = form.get_vertex_attribute(key, 'weight')
            dist += abs(z - targ)*weight

    return dist

def replicate_contraints(file, file_constraint):

    form = FormDiagram.from_json(file)
    form_ = FormDiagram.from_json(file_constraint)
    gkey_planar = {}

    for key_real in form.vertices():
        coord = form.vertex_coordinates(key_real)
        gkey_proj = geometric_key([coord[0],coord[1],0.0])
        gkey_planar[gkey_proj] = key_real

    for key in form_.vertices():
        target = form_.vertex[key].get('target', 0.0)
        lb = form_.vertex[key].get('lb', 0.0)
        ub = form_.vertex[key].get('ub', 0.0)
        if target < 10**(-4):
            target = 0.00
        gkey = geometric_key([form_.vertex_coordinates(key)[0],form_.vertex_coordinates(key)[1], 0.0])
        form.set_vertex_attribute(gkey_planar[gkey], 'target', target)
        form.set_vertex_attribute(gkey_planar[gkey], 'lb', lb)
        form.set_vertex_attribute(gkey_planar[gkey], 'ub', ub)

    return form

def interp_surf(form):

    x = []
    y = []
    s = []

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            x.append(vertex.get('x'))
            y.append(vertex.get('y'))
            s.append(vertex.get('target'))

    from scipy import interpolate

    surf = interpolate.interp2d(x, y, s, kind = 'linear')

    return surf

def null_edges(form, plot=False):

    null_edges = []
    all_edges = []

    for u, v in form.edges():
        if form.get_edge_attribute((u,v), 'is_external') == False and form.get_edge_attribute((u,v), 'is_edge') == True:
            activ = 0
            coord_u = form.vertex_coordinates(u)
            coord_v = form.vertex_coordinates(v)
            ux = round(coord_u[0],3)
            uy = round(coord_u[1],3)
            vx = round(coord_v[0],3)
            vy = round(coord_v[1],3)
            mid_x, mid_y, _ = form.edge_midpoint(u,v)
            if uy == vy and ((uy is not 10.0 and vy is not 10.0) or (uy is not 0.0 and vy is not 0.0)):
                if (mid_y > mid_x and mid_y < 10 - mid_x) or (mid_y < mid_x and mid_y > 10 - mid_x):
                    if uy == 5.0 and vy == 5.0 and ux > 0.01 and vx > 0.01 and ux < 9.99 and vx < 9.99: # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u,v))
                        activ += 1
            if ux == vx and ((ux is not 10.0 and vx is not 10.0) or (ux is not 0.0 and vx is not 0.0)):
                if (mid_y > mid_x and mid_y > 10 - mid_x) or (mid_y < mid_x and mid_y < 10 - mid_x):
                    if ux == 5.0 and vx == 5.0 and uy > 0.01 and vy > 0.01 and uy < 9.99 and vy < 9.99: # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u,v))
                        activ += 1
            if activ == 0:
                all_edges.append((u,v))

    if plot:
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges(all_edges)
        plotter.draw_edges(null_edges, color='ff0000')
        plotter.show()

    return null_edges

def set_height_constraint(form, zmax = 1.0):

    for key in form.vertices():
        form.set_vertex_attribute(key, 'ub', value = zmax)

    return form

def set_cross_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb=False, weights = False, thk = None, tol = 0.00, set_heights = False, update_loads = False, t = 10.0, b = 1.0):

    """ Set Cross-Vault heights.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    ub_lb : bool (optional)
        If True, the thickness will apply and the limits will be stored as attributes 'ub' and 'lb' on the form-diagram

    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface

    tol : float (optional)
        Approximates the equations avoiding negative square-roots.

    set_heights: bool
        If True, the nodes will have the heights 'z' updated to match the pointed arch shape.

    update_loads: bool
        If True, loads will be redistributed to match the target shape (Note: the magnitude will continue the  as in the input FD).

    t: float
        Negative lower-bound for the reactions position.

    Returns
    -------
    obj
        FormDiagram.

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4

    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx,ry)

    form_ = deepcopy(form)
    pzt = 0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
        yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
        hxd = (math.sqrt((rx)**2 - ((xd - x0)- rx)**2))
        hyd = (math.sqrt((ry)**2 - ((yd - y0)- ry)**2))
        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q1
            z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q3
            z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q2
            z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q4
            z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target', value = z)
        form_.set_vertex_attribute(key, 'z', value = z)
        pzt += form.get_vertex_attribute(key, 'pz')
        if weights:
            form.set_vertex_attribute(key,name='weight',value = form.get_vertex_attribute(key,'pz'))
        if set_heights:
            form.set_vertex_attribute(key, 'z', value = round(z,2))
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            form.set_vertex_attribute(key, 'b', value = [b,b])

    if update_loads:
        pz_3d = 0
        for key in form.vertices():
            pzi = form_.vertex_area(key)
            form.set_vertex_attribute(key, 'pz', value = pzi)
            pz_3d += pzi
        factor = pzt/pz_3d
        print('Planar load: {0:.2f} \ 3d Load: {1:.2f} \ Factor applied: {2:.3f}'.format(pzt,pz_3d,factor))
        for key in form.vertices():
            pz = factor * form.get_vertex_attribute(key, 'pz')
            form.set_vertex_attribute(key, 'pz', value = pz)


    if ub_lb:
        form = set_cross_vault_heights_ub(form, xy_span = xy_span, thk = thk)
        form = set_cross_vault_heights_lb(form, xy_span = xy_span, thk = thk, t = t)

    return form

def set_cross_vault_heights_ub(form, xy_span = [[0.0,10.0],[0.0,10.0]], thk = 0.5, tol = 0.000, set_heights = False):

    y1 = xy_span[1][1] + thk/2
    y0 = xy_span[1][0] - thk/2
    x1 = xy_span[0][1] + thk/2
    x0 = xy_span[0][0] - thk/2

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx,ry)

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
        yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
        hxd = (math.sqrt((rx)**2 - ((xd - x0)- rx)**2))
        hyd = (math.sqrt((ry)**2 - ((yd - y0)- ry)**2))
        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q1
            z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q3
            z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q2
            z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q4
            z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key, 'ub', value = z)
        if form.get_vertex_attribute(key, 'is_fixed'):
            form.attributes['tmax'] = z
        if set_heights:
            form.set_vertex_attribute(key,'z',value=z)

    return form

def set_cross_vault_heights_lb(form, xy_span = [[0.0,10.0],[0.0,10.0]], thk = 0.5, tol = 0.0, set_heights = False, t = 3.0):

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx,ry)

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if ((yi) > y1 and ((xi) > x1 or (xi) < x0 )) or ((yi) < y0 and ((xi) > x1 or (xi) < x0)):
            z = - 1* max(t,max(rx,ry)/4)
        else:
            if yi > y1:
                yi = y1
            elif yi < y0:
                yi = y0
            elif xi > x1:
                xi = x1
            elif xi < x0:
                xi = x0
            xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
            yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
            hxd = (_sqrt((rx)**2 - ((xd - x0) - rx)**2))
            hyd = (_sqrt((ry)**2 - ((yd - y0) - ry)**2))
            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q1
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q3
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q2
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q4
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            else:
                print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
                z = 0.0
        form.set_vertex_attribute(key,'lb',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=z)

    return form

def set_pavillion_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb = False, thk = 0.5, tol = 0.00, set_heights = False, b = 5.0, t = 2.0, update_loads = False, delete_boundary_edges = False):

    # Uodate this function to work on rectangular vaults

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2

    pzt = 0

    form_ = deepcopy(form)


    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q1
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q3
            z = math.sqrt((rx)**2 - (yi-rx)**2)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q2
            z = math.sqrt((ry)**2 - (xi-ry)**2)
        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q4
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key, 'target', value = z)
        form_.set_vertex_attribute(key, 'z', value = z)
        pzt += form.get_vertex_attribute(key, 'pz')
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            form.set_vertex_attribute(key, 'b', value = [b,b]) # Change to be [b,0] or [0,b]

    if update_loads:
        pz_3d = 0
        for key in form.vertices():
            pzi = form_.vertex_area(key)
            form.set_vertex_attribute(key, 'pz', value = pzi)
            pz_3d += pzi
        factor = pzt/pz_3d
        print('Planar load: {0:.2f} \ 3d Load: {1:.2f} \ Factor applied: {2:.3f}'.format(pzt,pz_3d,factor))
        for key in form.vertices():
            pz = factor * form.get_vertex_attribute(key, 'pz')
            form.set_vertex_attribute(key, 'pz', value = pz)

    if ub_lb:
        form = set_pavillion_vault_heights_ub(form, xy_span = xy_span, thk = thk)
        form = set_pavillion_vault_heights_lb(form, xy_span = xy_span, thk = thk, t = t)

    return form

def set_pavillion_vault_heights_ub(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb = False, thk = 0.5, tol = 0.00, set_heights = False):

    y1 = xy_span[1][1] + thk/2
    y0 = xy_span[1][0] - thk/2
    x1 = xy_span[0][1] + thk/2
    x0 = xy_span[0][0] - thk/2

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] + thk)/2.0
    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
            z = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
            z = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
            z = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
            z = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'ub',value=z)
        if form.get_vertex_attribute(key, 'is_fixed'):
            form.attributes['tmax'] = z
        # form.set_vertex_attribute(key,'z',value=z)

    return form

def set_pavillion_vault_heights_lb(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb = False, thk = 0.5, tol = 0.00, set_heights = False, t = 3.0):

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] - thk)/2.0


    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if (yi) > y1 or (xi) > x1 or (xi) < x0 or (yi) < y0:
            z = -1.0 * t
        else:
            if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
                z = math.sqrt(abs((ry)**2 - ((yi - y0)-ry)**2))
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
                z = math.sqrt(abs((ry)**2 - ((yi - y0)-ry)**2))
            elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
                z = math.sqrt(abs((rx)**2 - ((xi - x0)-rx)**2))
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
                z = math.sqrt(abs((rx)**2 - ((xi - x0)-rx)**2))
            else:
                print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
                z = 0.0
        form.set_vertex_attribute(key,'lb',value=z)
        # form.set_vertex_attribute(key,'z',value=z)

    return form


def set_oct_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], thickness = None, tol = 0.01):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi < y1/x1*xi + tol and yi < y1 - x1 + tol: #Q1
            z = math.sqrt= ((rx)**2 - (xi-rx)**2)
        elif yi > y1/x1*xi - tol and yi > y1 - x1 - tol: #Q3
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi < y1/x1*xi + tol and yi > y1 - x1 - tol: #Q2
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        elif yi > y1/x1*xi - tol and yi < y1 - x1 + tol: #Q4
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)

    return form

def set_dome_heights(form, center = [0.0,0.0], radius = 10.0, thck = 0.30, set_heights = False, t = 5.0, update_loads = True):

    x0 = center[0]
    y0 = center[1]
    ri = radius - thck/2
    re = radius + thck/2
    loads = deepcopy(form)
    pzt = 0.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        zt2 = radius**2 - (xi - x0)**2 - (yi - y0)**2
        zi2 = ri**2 - (xi - x0)**2 - (yi - y0)**2
        ze2 = re**2 - (xi - x0)**2 - (yi - y0)**2
        # z = math.sqrt(zt2)
        ze = math.sqrt(ze2)
        if zi2 < 0.0:
            zi = - t
        else:
            zi = math.sqrt(zi2)
        if -0.001 <= zt2 <= 0.0:
            zt2 = 0.0
        try:
            z = math.sqrt(zt2)
        except:
            print(xi,yi)
            z=0

        loads.set_vertex_attribute(key, 'z', value=z)
        form.set_vertex_attribute(key, 'target', value=z)
        form.set_vertex_attribute(key, 'lb', value=zi)
        form.set_vertex_attribute(key, 'ub', value=ze)
        pzt += form.get_vertex_attribute(key, 'pz')

    pz3d = 0
    for key in form.vertices():
        pz = loads.vertex_area(key)
        form.set_vertex_attribute(key,'pz',value=pz)
        pz3d += pz

    scale = pzt/pz3d
    print('3D Load: {0:.2f}, Planar Load: {1:.2f} / Scale: {2:.2f}'.format(pz3d, pzt, scale))

    if update_loads:
        pzt = 0.0
        for key in form.vertices():
            form.vertex[key]['pz'] *= scale
            pzt += form.vertex[key]['pz']

    print('After Scaling Load: {0:.2f}'.format(pzt))

    for key in form.vertices_where({'is_fixed': True}):
        x, y, _ = form.vertex_coordinates(key)
        theta = math.atan2((y - y0), (x - x0))
        x_ = thck/2*math.cos(theta)
        y_ = thck/2*math.sin(theta)
        form.set_vertex_attribute(key, 'b', [x_, y_])

    return form

def _find_r_given_h_l(h,l):

    r = h**2/l + l/4

    return r

def circle_3points_xy(p1,p2,p3):

    x1 = p1[0]
    z1 = p1[1]
    x2 = p2[0]
    z2 = p2[1]
    x3 = p3[0]
    z3 = p3[1]

    x12 = x1 - x2
    x13 = x1 - x3
    z12 = z1 - z2
    z13 = z1 - z3
    z31 = z3 - z1
    z21 = z2 - z1
    x31 = x3 - x1
    x21 = x2 - x1

    sx13 = x1**2 - x3**2
    sz13 = z1**2 - z3**2
    sx21 = x2**2 - x1**2
    sz21 = z2**2 - z1**2

    f = ((sx13) * (x12) + (sz13) * (x12) + (sx21) * (x13) + (sz21) * (x13)) / (2 * ((z31) * (x12) - (z21) * (x13)))
    g = ((sx13) * (z12) + (sz13) * (z12) + (sx21) * (z13) + (sz21) * (z13)) / (2 * ((x31) * (z12) - (x21) * (z13)))
    c = - x1 ** 2 -  z1 ** 2 - 2 * g * x1 - 2 * f * z1
    h = - g
    k = - f
    r2 = h * h + k * k - c
    r = math.sqrt(r2)

    print('h: ',h,'k: ',k,'r: ',r)

    return h, k, r

def _sqrt(x):
    try:
        sqrt_x = math.sqrt(x)
    except:
        if x > -10e4:
            sqrt_x = math.sqrt(abs(x))
        else:
            sqrt_x = 0.0
            print('Problems to sqrt: ',x)
    return sqrt_x

def set_pointed_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], hc=8.0, he=None, hm=None,  ub_lb=False, thk = None, tol = 0.00, set_heights = False):

    """ Set pointed-vault heights.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    hc: float
        Height of the central part of the pointed vault.

    he: list (optional)
        Height of the opening mid-span for each of the quadrants (see Notes).

    hm: list (optional)
        Height of each quadrant mid-span (see Notes).

    ub_lb : bool (optional)
        If True, the thickness will apply and the limits will be stored as attributes 'ub' and 'lb' on the form-diagram

    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface

    tol : float (optional)
        Approximates the equations avoiding negative square-roots.

    set_heights: bool
        If True, the nodes will have the heights 'z' updated to match the pointed arch shape.

    Returns
    -------
    obj
        FormDiagram.

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4

    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    lx = x1 - x0
    ly = y1 - y0

    if he and hm is None:
        h1, k1, r1 = circle_3points_xy([x0,he[1]],[(x1+x0)/2,hc],[x1,he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0,he[3]],[(y1+y0)/2,hc],[y1,he[2]])
        h4, k4, r4 = h3, k3, r3
    elif hm and he:
        h1, k1, r1 = circle_3points_xy([(x1+x0)/2,hc],[3*(x1+x0)/4,hm[0]],[x1,he[0]])
        h2, k2, r2 = circle_3points_xy([(x1+x0)/2,hc],[1*(x1+x0)/4,hm[1]],[x0,he[1]])
        h3, k3, r3 = circle_3points_xy([(y1+y0)/2,hc],[3*(y1+y0)/4,hm[2]],[y1,he[2]])
        h4, k4, r4 = circle_3points_xy([(y1+y0)/2,hc],[1*(y1+y0)/4,hm[3]],[y0,he[3]])

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)

        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q1
            # Equation (xi - hx) ** 2 + (hi - kx) ** 2 = rx **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,ly) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q3
            # Equation (xi - hy) ** 2 + (hi - ky) ** 2 = ry **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k3 + math.sqrt(r3 ** 2 - (yi - h3) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,lx) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q2
            if he:
                hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,ly)
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q4
            if he:
                hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,lx)
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
        form.set_vertex_attribute(key,'target',value=zi)

        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(zi,2))

    if ub_lb and thk:
        form = set_pointed_vault_heights_ub(form, xy_span = xy_span, hc = hc, he = he, hm = hm, thk = thk)
        form = set_pointed_vault_heights_lb(form, xy_span = xy_span, hc = hc, he = he, hm = hm, thk = thk)

    return form

def set_pointed_vault_heights_ub(form, xy_span = [[0.0,10.0],[0.0,10.0]], hc=8.0, he=None, hm=None, thk = None, tol = 0.00, set_heights = False):

    y1_init = xy_span[1][1]
    y0_init = xy_span[1][0]
    x1_init = xy_span[0][1]
    x0_init = xy_span[0][0]

    y1 = xy_span[1][1] + thk / 2
    y0 = xy_span[1][0] - thk / 2
    x1 = xy_span[0][1] + thk / 2
    x0 = xy_span[0][0] - thk / 2

    lx = x1 - x0
    ly = y1 - y0

    hc += thk/2
    if he:
        for i in range(len(he)):
            he[i] += thk/2
    if hm:
        for i in range(len(he)):
            hm[i] += thk/2

    if he and hm is None:
        h1, k1, r1 = circle_3points_xy([x0_init,he[1]],[(x1_init+x0_init)/2,hc],[x1_init,he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0_init,he[3]],[(y1_init+y0_init)/2,hc],[y1_init,he[2]])
        h4, k4, r4 = h3, k3, r3
    elif hm and he:
        h1, k1, r1 = circle_3points_xy([(x1_init+x0_init)/2,hc],[3*(x1_init+x0_init)/4,hm[0]],[x1_init,he[0]])
        h2, k2, r2 = circle_3points_xy([(x1_init+x0_init)/2,hc],[1*(x1_init+x0_init)/4,hm[1]],[x0_init,he[1]])
        h3, k3, r3 = circle_3points_xy([(y1_init+y0_init)/2,hc],[3*(y1_init+y0_init)/4,hm[2]],[y1_init,he[2]])
        h4, k4, r4 = circle_3points_xy([(y1_init+y0_init)/2,hc],[1*(y1_init+y0_init)/4,hm[3]],[y0_init,he[3]])

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)

        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q1
            # Equation (xi - hx) ** 2 + (hi - kx) ** 2 = rx **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,ly) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q3
            # Equation (xi - hy) ** 2 + (hi - ky) ** 2 = ry **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k3+ math.sqrt(r3 ** 2 - (yi - h3) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,lx) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q2
            if he:
                hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,ly)
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q4
            if he:
                hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi,lx)
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
        form.set_vertex_attribute(key,'target',value=zi)

        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(zi,2))

    return form

def set_pointed_vault_heights_lb(form, xy_span = [[0.0,10.0],[0.0,10.0]], hc=8.0, he=None, hm=None, thk = None, tol = 0.00, set_heights = False):

    y1_init = xy_span[1][1]
    y0_init = xy_span[1][0]
    x1_init = xy_span[0][1]
    x0_init = xy_span[0][0]

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    lx = x1 - x0
    ly = y1 - y0

    hc -= thk/2
    if he:
        for i in range(len(he)):
            he[i] -= thk/2
    if hm:
        for i in range(len(he)):
            hm[i] -= thk/2

    if hm and he:
        h1, k1, r1 = circle_3points_xy([(x1_init+x0_init)/2,hc],[3*(x1_init+x0_init)/4,hm[0]],[x1_init,he[0]])
        h2, k2, r2 = circle_3points_xy([(x1_init+x0_init)/2,hc],[1*(x1_init+x0_init)/4,hm[1]],[x0_init,he[1]])
        h3, k3, r3 = circle_3points_xy([(y1_init+y0_init)/2,hc],[3*(y1_init+y0_init)/4,hm[2]],[y1_init,he[2]])
        h4, k4, r4 = circle_3points_xy([(y1_init+y0_init)/2,hc],[1*(y1_init+y0_init)/4,hm[3]],[y0_init,he[3]])
    elif he and hm is None:
        h1, k1, r1 = circle_3points_xy([x0_init,he[1]],[(x1_init+x0_init)/2,hc],[x1_init,he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0_init,he[3]],[(y1_init+y0_init)/2,hc],[y1_init,he[2]])
        h4, k4, r4 = h3, k3, r3

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)

        if (yi) > y1 or (xi) > x1 or (xi) < x0 or (yi) < y0:
            zi = 0.0
        else:

            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q1
                # Equation (xi - hx) ** 2 + (hi - kx) ** 2 = rx **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
                if he:
                    hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi,ly) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
                if yi <= (y1 + y0)/2:
                    zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q3
                # Equation (xi - hy) ** 2 + (hi - ky) ** 2 = ry **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
                if he:
                    hi = k3+ math.sqrt(r3 ** 2 - (yi - h3) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi,lx) # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
                if xi <= (x0 + x1)/2:
                    zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: # Q2
                if he:
                    hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi,ly)
                if yi <= (y1 + y0)/2:
                    zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: # Q4
                if he:
                    hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi,lx)
                if xi <= (x0 + x1)/2:
                    zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

            else:
                print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))

        form.set_vertex_attribute(key,'target',value=zi)

        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(zi,2))

    return form

def circular_heights(form , x0 = None, xf = None, thk = 0.5, t = 20, density = 20, ly = 1.0, overwrite_weight=False):
    """ Create constraints for a circular arch.

    Parameters
    ----------
    form : FormDiagram
        linear form diagram to apply he constraints.

    x0: float (optional)
        Beginning of the linear form diagram.

    xf: float (optional)
        End of the linear form diagram.

    thk : float
        Homothetical (extruded to both sides) thickness of arch.

    t: float
        Parameter to constraint from below the heights of the restraint and/or the nodes that have no vertical projection on arch's intrados.

    density: float
        Density of the masonry, 20 kN/m3 is taken by default, if not desired put one.

    ly : float
        Out of plane dimension of the arch. Default = 1.0

    Returns
    -------
    obj
        FormDiagram.

    """


    if x0 == None or xf == None:
        x = []
        for key in form.vertices():
            x.append(form.vertex_coordinates(key)[0])
        x0 = min(x)
        xf = max(x)

    xc = (xf+x0)/2
    r = xf - xc
    ri = r - thk/2
    re = r + thk/2
    form.attributes['Re'] = re
    form.attributes['Ri'] = ri
    if overwrite_weight is False:
        arch_area = (math.pi*re**2 - math.pi*ri**2)/2
        total_weight = arch_area * density * ly
        weight_unit = total_weight / (form.number_of_vertices() - 1)
    else:
        total_weight = overwrite_weight
        weight_unit = total_weight / (form.number_of_vertices() - 1)

    print('SpanMid: {0:.2} m / SpanInt: {1:.2} m / SpanExt: {2:.2} m / Thickness: {3:.4} m / Ratio t/Ri: {4:.4} m / Ratio t/R: {5:.4} m / Weight: {6:.3}'.format(2*r, 2*ri, 2*re, thk, (thk/ri), (thk/r), total_weight))
    pzt = 0

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        zt = math.sqrt(r**2 - (x-xc)**2)
        ze = math.sqrt(re**2 - (x-xc)**2)
        form.set_vertex_attribute(key,'target',value=zt)
        zi2 = ri**2 - (x-xc)**2
        if zi2 < 0:
            zi = 0 - t
        else:
            zi = math.sqrt(zi2)
        if form.get_vertex_attribute(key,'is_fixed') == True:
            form.set_vertex_attribute(key,'lb',value=-t)
            form.set_vertex_attribute(key,'ub',value=ze)
        else:
            form.set_vertex_attribute(key,'lb',value=zi)
            form.set_vertex_attribute(key,'ub',value=ze)
        # form.set_vertex_attribute(key,'z',value=ze)
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            form.set_vertex_attribute(key, 'b', value = [thk/2,0.0])
        if x == x0:
            form.attributes['tmax'] = ze
        if form.get_vertex_attribute(key, 'is_fixed') is True:
            pi = weight_unit / 2
            form.set_vertex_attribute(key, 'pz', pi)
            pzt += pi
        else:
            pi = weight_unit
            form.set_vertex_attribute(key, 'pz', pi)
            pzt += pi
    print('Total weight applied: {0:.3}'.format(pzt))
    form.attributes['selfweight'] = pzt

    return form

def create_cracks(form , dx =[[0.50, 0.55]], dy = [[-0.1, 0.1]], type = ['top'], view = False):

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

def circular_joints(form , x0 = None, xf = None, blocks = 18, thk=0.5, t=0.0, tol = 1e-3):

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
    print('SpanMid: {0:.2} m / SpanInt: {1:.2} m / SpanExt: {2:.2} m / Thickness: {3:.4} m / Ratio t/Ri: {4:.4} m / Ratio t/R: {5:.4} m / Number of Blocks: {6}'.format(2*r, 2*ri, 2*re, thk, (thk/ri), (thk/r),blocks))

    njoints = blocks+1
    joints = {}
    for j in range(njoints):
        theta = j/blocks*math.pi
        xi = xc + ri * math.cos(theta) # takeout
        zi = ri * math.sin(theta) + tol
        xe = xc + re * math.cos(theta)
        ze = re * math.sin(theta) + tol # take out
        xmax = max(xi,xe)
        xmin = min(xi,xe)
        possible_edges = []
        for (u,v) in form.edges():
            xu, xv = form.vertex_coordinates(u)[0], form.vertex_coordinates(v)[0]
            if max(xu,xv) >= xmin and min(xu,xv) <= xmax:
                possible_edges.append(tuple(sorted([k_i[u],k_i[v]])))
                if form.get_vertex_attribute(u, 'is_fixed') == True:
                    possible_edges.append(tuple(sorted([-k_i[u],k_i[u]])))
                if form.get_vertex_attribute(v, 'is_fixed') == True:
                    possible_edges.append(tuple(sorted([-k_i[v],k_i[v]])))
        joints[j] = [[xi,y,zi],[xe,y,ze],set(possible_edges)]
        print(joints[j])
    form.attributes['joints'] = joints

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        zt = math.sqrt(r**2 - (x-xc)**2)
        ze = math.sqrt(re**2 - (x-xc)**2) - t
        form.set_vertex_attribute(key,'target',value=zt)
        zi2 = ri**2 - (x-xc)**2
        if zi2 < 0:
            zi = 0 - t
        else:
            zi = math.sqrt(zi2) - t
        if form.get_vertex_attribute(key,'is_fixed') == True:
            form.set_vertex_attribute(key,'lb',value=None)
            form.set_vertex_attribute(key,'ub',value=None)
        else:
            form.set_vertex_attribute(key,'lb',value=zi)
            form.set_vertex_attribute(key,'ub',value=ze)
        # form.set_vertex_attribute(key,'z',value=ze)
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            form.set_vertex_attribute(key, 'b', value = [thk/2,0.0])
        if x == x0:
            form.attributes['tmax'] = ze

    return form

def rollers_on_openings(form, xy_span = [[0.0,10.0],[0.0,10.0]], max_f = 5.0, constraint_directions = 'all'):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    bndr = form.vertices_on_boundary()

    for key in bndr:
        if form.get_vertex_attribute(key, 'is_fixed') == False:
            x, y, _ = form.vertex_coordinates(key)
            if x == x1 and (constraint_directions in ['all', 'x']):
                form.set_vertex_attribute(key, 'rol_x', True)
                form.set_vertex_attribute(key, 'max_rx', max_f)
            if x == x0 and (constraint_directions in ['all', 'x']):
                form.set_vertex_attribute(key, 'rol_x', True)
                form.set_vertex_attribute(key, 'max_rx', max_f)
            if y == y1 and (constraint_directions in ['all', 'y']):
                form.set_vertex_attribute(key, 'rol_y', True)
                form.set_vertex_attribute(key, 'max_ry', max_f)
            if y == y0 and (constraint_directions in ['all', 'y']):
                form.set_vertex_attribute(key, 'rol_y', True)
                form.set_vertex_attribute(key, 'max_ry', max_f)

    return form
