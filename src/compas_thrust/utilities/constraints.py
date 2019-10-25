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
from compas_thrust.algorithms.equilibrium import z_from_form

from compas_plotters import MeshPlotter

import math

from compas_thrust.plotters.plotters import plot_form

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
    'circular_heights'
]


def check_constraints(form, show=False, lb_show=False, ub_show=False):

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
            if z < lb:
                outside[key] = lb - z
                penalty += (abs(outside[key])+4)**(4)
        if vertex.get('ub', None):
            ub = vertex['ub']
            ubs[key] = ub
            if z > ub:
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

def set_cross_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb=False, weights = False, thk = None, tol = 0.00, set_heights = False):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi <= y1/x1 * xi + tol and yi <= y1 - xi + tol: #Q1
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi >= y1/x1 * xi - tol and yi >= y1 - xi - tol: #Q3
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi <= y1/x1 * xi + tol and yi >= y1 - xi - tol: #Q2
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        elif yi >= y1/x1 * xi - tol and yi <= y1 - xi + tol: #Q4
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)
        if weights:
            form.set_vertex_attribute(key,name='weight',value = form.get_vertex_attribute(key,'pz'))
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    if ub_lb:
        form = set_cross_vault_heights_ub(form, xy_span = xy_span, thk = thk)
        form = set_cross_vault_heights_lb(form, xy_span = xy_span, thk = thk)

    return form

def set_cross_vault_heights_ub(form, xy_span = [[0.0,10.0],[0.0,10.0]], thk = 0.5, tol = 0.000):

    y1 = xy_span[1][1] + thk/2
    y0 = xy_span[1][0] - thk/2
    x1 = xy_span[0][1] + thk/2
    x0 = xy_span[0][0] - thk/2

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] + thk)/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
            z = math.sqrt(((rx)**2 - ((xi - x0)-rx)**2))
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
            z = math.sqrt(((rx)**2 - ((xi - x0)-rx)**2))
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
            z = math.sqrt(((ry)**2 - ((yi - y0)-ry)**2))
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
            z = math.sqrt(((ry)**2 - ((yi - y0)-ry)**2))
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key, 'ub', value = z)

    return form

def set_cross_vault_heights_lb(form, xy_span = [[0.0,10.0],[0.0,10.0]], thk = 0.5, tol = 0.0):

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] - thk)/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if ((yi) > y1 and ((xi) > x1 or (xi) < x0 )) or ((yi) < y0 and ((xi) > x1 or (xi) < x0)):
            z = 0.0
        else:
            if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
                z = math.sqrt(((rx)**2 - ((xi - x0)-rx)**2))
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
                z = math.sqrt(((rx)**2 - ((xi - x0)-rx)**2))
            elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
                z = math.sqrt(((ry)**2 - ((yi - y0)-ry)**2))
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
                z = math.sqrt(((ry)**2 - ((yi - y0)-ry)**2))
            else:
                print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
                z = 0.0
        form.set_vertex_attribute(key,'lb',value=z)

    return form

def set_pavillion_vault_heights(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb = False, thk = 0.5, tol = 0.00, set_heights = False):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi <= y1/x1 * xi + tol and yi <= y1 - xi + tol: #Q1
            try:
                z = math.sqrt((ry)**2 - (yi-ry)**2)
            except:
                z = 0
        elif yi >= y1/x1 * xi - tol and yi >= y1 - xi - tol: #Q3
            try:
                z = math.sqrt((ry)**2 - (yi-ry)**2)
            except:
                z = 0
        elif yi <= y1/x1 * xi + tol and yi >= y1 - xi - tol: #Q2
            try:
                z = math.sqrt((rx)**2 - (xi-rx)**2)
            except:
                z = 0
        elif yi >= y1/x1 * xi - tol and yi <= y1 - xi + tol: #Q4
            try:
                z = math.sqrt((rx)**2 - (xi-rx)**2)
            except:
                z = 0
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    if ub_lb:
        form = set_pavillion_vault_heights_ub(form, xy_span = xy_span, thk = thk)
        form = set_pavillion_vault_heights_lb(form, xy_span = xy_span, thk = thk)

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
        # form.set_vertex_attribute(key,'z',value=z)

    return form

def set_pavillion_vault_heights_lb(form, xy_span = [[0.0,10.0],[0.0,10.0]], ub_lb = False, thk = 0.5, tol = 0.00, set_heights = False):

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] - thk)/2.0
    

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if (yi) > y1 or (xi) > x1 or (xi) < x0 or (yi) < y0:
            z = 0.0
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
        form.set_vertex_attribute(key,'z',value=z)

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

def set_dome_heights(form, center = [0.0,0.0], radius = 10.0, thickness = None, tol = 0.00, set_heights = False):

    x0 = center[0]
    y0 = center[1]

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        z2 = + radius**2 - (xi - x0)**2 - (yi - y0)**2
        if -0.01 <= z2 <= 0.0:
            z2 = 0.0
        try: 
            z = math.sqrt(z2)
        except:
            print(xi,yi)
            z=0
        form.set_vertex_attribute(key,'target',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    return form


def set_dome_heights(form, center = [0.0,0.0], radius = 10.0, thickness = None, tol = 0.00, set_heights = False):

    x0 = center[0]
    y0 = center[1]

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        z2 = + radius**2 - (xi - x0)**2 - (yi - y0)**2
        if -0.01 <= z2 <= 0.0:
            z2 = 0.0
        try: 
            z = math.sqrt(z2)
        except:
            print(xi,yi)
            z=0
        form.set_vertex_attribute(key,'target',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    return form

def circular_heights(form , x0 = None, xf = None, thk=0.5, t=0.0):

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
    print('SpanMid: {0:.2} m / SpanInt: {1:.2} m / Thickness: {2:.2} m / Ratio t/Ri: {3:.2} m'.format(2*xc, 2*ri, thk, (thk/ri)))

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
        if x == x0:
            form.attributes['tmax'] = ze
        

    return form