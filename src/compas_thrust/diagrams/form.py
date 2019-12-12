from compas.utilities import geometric_key
from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax
from compas.numerical import connectivity_matrix

from compas_thrust.algorithms.equilibrium import z_from_form
from compas_thrust.algorithms.scale import scale_form
from compas_thrust.algorithms.scale import evaluate_scale
from compas_thrust.algorithms.equilibrium import z_update


from compas_thrust.plotters.plotters import plot_form
from compas_thrust.plotters.plotters import plot_force

from random import shuffle
from copy import copy
from compas_tna.utilities import parallelise_sparse

from numpy.random import rand
from numpy import array
from numpy import float64
from numpy import zeros
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum

import math


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    '_form',
    'adapt_tna',
    'adapt_objective',
    'evaluate_a',
    'remove_feet',
    'energy',
    'loadpath',
    'delete_boundary_edges',
    'overview_forces',
    'create_arch',
    'create_dome_form',
]

def _form(form, keep_q=False):

    """ s the FormDiagram by shuffling the edges.

    Parameters
    ----------
    form : obj
        Original FormDiagram.

    Returns
    -------
    obj
        Shuffled FormDiagram.

    """

    # Edges

    edges = [form.edge_coordinates(u, v) for u, v in form.edges()]
    edges = [[sp[:2] + [0], ep[:2] + [0]] for sp, ep in edges]
    qs = {geometric_key(form.edge_midpoint(u, v)[:2] + [0]) : form.get_edge_attribute((u,v), 'q') for u, v in form.edges()}
    shuffle(edges)

    form_ = FormDiagram.from_lines(edges, delete_boundary_face=False)
    form_.update_default_edge_attributes({'is_symmetry': False})
    sym = [geometric_key(form.edge_midpoint(u, v)[:2] + [0])for u, v in form.edges_where({'is_symmetry': True})]
    for u, v in form_.edges():
        if geometric_key(form_.edge_midpoint(u, v)) in sym:
            form_.set_edge_attribute((u, v), 'is_symmetry', True)
        if keep_q:
            form_.set_edge_attribute((u, v), 'q', qs[geometric_key(form_.edge_midpoint(u, v)[:2] + [0])])

    # Vertices

    gkey_key = form_.gkey_key()
    for key, vertex in form.vertex.items():
        gkey = geometric_key(form.vertex_coordinates(key)[:2] + [0])
        form_.vertex[gkey_key[gkey]] = vertex

    form_.attributes['indset'] = []

    return form_

def add_feet(form, delete_face = False, plot = False):

    if delete_face:
        form.delete_face(0)
    corners = list(form.vertices_where({'is_fixed': True}))
    form.set_vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
    form.update_boundaries(feet=2)
    if plot:
        form.plot()

    return form

def adapt_tna(form, zmax = 5.0, method = 'nodal', plot = False, delete_face = False, alpha = 100.0, kmax = 100, display = False):

    form = add_feet(form, delete_face = delete_face, plot = plot)
    force  = ForceDiagram.from_formdiagram(form)

    if method == 'nodal':
        horizontal_nodal(form, force, alpha = alpha, kmax = kmax, display = False)
    else:
        horizontal(form, force, alpha=alpha, kmax = kmax, display = False)

    # Vertical Equilibrium with no updated loads

    form0 = copy(form)
    form_ = scale_form(form0, 1.0)
    z = [form_.get_vertex_attribute(key, 'z') for key in form_.vertices()]
    z_ = max(z)
    scale = (z_ / zmax)
    form = scale_form(form0, scale)
    
    if plot:
        plot_form(form).show()
        plot_force(force, form).show()
        force.plot()

    return form

def adapt_objective(form, zrange = [3.0,8.0], objective = 'loadpath', method = 'nodal', discr = 100, plot = False, 
                    delete_face = False, alpha = 100.0, kmax = 100, display = False, amax = 2.0, rmax = 0.01):

    form = add_feet(form, delete_face = delete_face, plot = plot)
    force  = ForceDiagram.from_formdiagram(form)

    plot_force(force, form).show()

    if method == 'nodal':
        horizontal_nodal(form, force, alpha = alpha, kmax = kmax, display = False)
    else:
        horizontal(form, force, alpha=alpha, kmax = kmax, display = False)

    plot_force(force, form).show()

    a = evaluate_a(form, plot = plot)
    if a > amax:
        print('High Angle deviations!')

    # Vertical Equilibrium with no updated loads

    if objective == 'loadpath':
        f = loadpath
    if objective == 'target':
        f = energy
    
    scale = []
    form0 = scale_form(form, 1.0)
    z = [form0.get_vertex_attribute(key, 'z') for key in form0.vertices()]
    z_ = max(z)
    scale.append(z_ / zrange[1])
    scale.append(z_ / zrange[0])
    print(scale)

    best_scl = evaluate_scale(form0,f,scale,n=discr, plot = plot)
    
    print('Best Scale is: {0}'.format(best_scl))
    form = scale_form(form0, best_scl)

    r = residual(form, plot = plot)
    if r > rmax:
        print('High residual forces!')
    
    if plot:
        plot_form(form).show()
        plot_force(force, form).show()

    return form

def remove_feet(form, plot = False, openings = None, rmax = 0.01): # Flatten Diagram

    lines = []
    qs = {}

    if plot:
        plot_form(form).show()

    for u, v in form.edges_where({'is_edge': True, 'is_external': False}):
        s = form.vertex_coordinates(u)
        e = form.vertex_coordinates(v)
        lines.append([s,e])
        qs[geometric_key(form.edge_midpoint(u,v))] = form.get_edge_attribute((u,v), 'q')

    fixed = [geometric_key(form.vertex_coordinates(key)) for key in form.vertices_where({'is_anchor': True })]
    zs = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.vertex_coordinates(key)[2] for key in form.vertices_where({'is_external': False })}
    pz = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'pz') for key in form.vertices()}
    target = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'target') for key in form.vertices()}
    lb = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'lb') for key in form.vertices()}
    ub = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'ub') for key in form.vertices()}
    
    form_ = FormDiagram.from_lines(lines)
    form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False, 'is_edge': True})
    form_.update_default_vertex_attributes({'is_roller': False})

    if openings:
        for key in form.faces():
            if form.face_area(key) > openings - 1.0 and form.face_area(key) < openings + 1.0:
                form.delete_face(key)
                print('Deleted area of face {0}'.format(key))
                break
    gkey_key = form_.gkey_key()

    for pt in fixed:
        form_.set_vertex_attribute(gkey_key[pt], name = 'is_fixed', value = True)

    for key, attr in form_.vertices(True):
        pzi = pz[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
        zi = zs[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
        ti = target[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
        ub_i = ub[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
        lb_i = lb[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
        attr['pz'] = pzi
        attr['z'] = zi
        attr['target'] = ti
        attr['lb'] = lb_i
        attr['ub'] = ub_i

    for u, v in form_.edges():
        qi = qs[geometric_key(form_.edge_midpoint(u,v))]
        form_.set_edge_attribute((u,v), name = 'q', value = qi)
        print(qi)

    # form_ = z_from_form(form_) # This moves also x,y a bit...
    form_ = z_update(form_)

    r = residual(form, plot = plot)
    if r > rmax:
        print('High residual forces!')

    if plot:
        plot_form(form_).show()

    return form_

def evaluate_a(form, plot=True):

    a_total = 0
    a_max = 0
    for u, v, attr in form.edges_where({'is_edge': True}, True):
        a = attr['a']
        a_total += a
        l = form.edge_length(u,v)
        a = a*l
        if a > a_max:
            a_max = a
    if plot is True:
        print('Angle Deviation  Max: {0}'.format(a_max))
        print('Angle Deviation  Total: {0}'.format(a_total))

    return a_max

def residual(form, plot = False):

    # Mapping

    k_i  = form.key_index()

    # Vertices and edges

    n     = form.number_of_vertices()
    fixed = [k_i[key] for key in form.fixed()]
    rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
    free  = list(set(range(n)) - set(fixed) - set(rol))
    

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)

    xy = xyz[:, :2]
    px = px[free]
    py = py[free]
    pz = pz[free]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cit = Ci.transpose()
    uvw = C.dot(xyz)
    U   = uvw[:, 0]
    V   = uvw[:, 1]
    q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

    # Horizontal checks

    Rx = Cit.dot(U * q.ravel()) - px.ravel()
    Ry = Cit.dot(V * q.ravel()) - py.ravel()
    R = sqrt(Rx**2 + Ry**2)
    Rmax  = max(R)
    Rs = sum(R)

    if plot:
        print('Residual Total: {0}'.format(Rs))
        print('Residual Max: {0}'.format(Rmax))
    
    return Rmax 

def energy(form):

    f = 0
    for key, vertex in form.vertex.items():
            if vertex.get('is_external') == False:
                z  = vertex.get('z')
                s  = vertex.get('target')
                w = vertex.get('weight', 1.0)
                f += w * (z - s)**2

    return f

def loadpath(form):

    lp = 0
    for u, v in form.edges_where({'is_external': False}):
            if form.get_edge_attribute((u,v),'is_edge') is True and form.get_edge_attribute((u,v),'is_symmetry') is False:
                qi = form.get_edge_attribute((u,v),'q')
                li = form.edge_length(u,v)
                lp += qi*li**2

    return lp

def delete_boundary_edges(form):

    for u,v in form.edges_on_boundary():
        form.set_edge_attribute((u,v), 'is_edge', False)

    return form

def overview_forces(form):

    f = []
    q = []
    z = []
    pz = 0

    lp=0

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True and form.get_edge_attribute((u,v),'is_symmetry') is False:
            qi = form.get_edge_attribute((u,v),'q')
            li = form.edge_length(u,v)
            lp += qi*li**2
            q.append(qi)
            f.append(qi*li)

    print('='*20)
    print('Overview on forces:')

    print('q: {0:.3f} : {1:.3f}'.format(float(min(q)), float(max(q))))
    print('f: {0:.3f} : {1:.3f}'.format(float(min(f)), float(max(f))))
    for key in form.vertices():
        z.append(form.get_vertex_attribute(key,'z'))
        pz += form.get_vertex_attribute(key,'pz')
    print('z: {0:.3f} : {1:.3f}'.format(float(min(z)), float(max(z))))
    print('pz: {0:.3f}'.format(pz))
    print('lp: {0:.3f}'.format(lp))
    form.attributes['loadpath'] = lp

    return

def create_arch(D = 2.00, x0 = 0.0, total_nodes = 100, total_self_weight = 20.0, lambda_hor = None, f_concentrated = None ):

    r = D/2
    xc = x0 + r
    lines = []
    total_edges = total_nodes - 1
    gkey_fix = []

    for i in range(total_edges):
        xi = xc - r*math.cos(i/total_edges*math.pi)
        xf = xc - r*math.cos((i+1)/total_edges*math.pi)
        lines.append([[xi,0.0,0.0],[xf,0.0,0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi,0.0,0.0], precision=6))
        elif i == total_edges - 1:
            gkey_fix.append(geometric_key([xf,0.0,0.0], precision=6))

    form = FormDiagram.from_lines(lines, delete_boundary_face=False)
    gkey_key = form.gkey_key(precision=6)

    form.update_default_vertex_attributes({'is_roller': False})
    form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form.attributes['loadpath'] = 0
    form.attributes['indset'] = []
    form.set_vertices_attribute('pz', total_self_weight/(total_nodes))
    if lambda_hor:
        form.set_vertices_attribute('px', lambda_hor*total_self_weight/(total_nodes))
    if  f_concentrated:
        lbd, key, direc = f_concentrated
        form.set_vertex_attribute(key,direc,value=lbd*form.get_vertex_attribute(key,'pz'))   
    form.set_vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.set_vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form


def create_cross_form(xy_span = [[0.0,10.0],[0.0,10.0]], division = 10, fix = 'corners', rollers = False ,pz = True, px = None, py = None):

    """ Create a cross form-diagram.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    division: int
        Set the density of the grid in x and y directions.

    fix : string
        Option to select the constrained nodes: 'corners', 'all' are accepted.

    rollers: bool
        Set the non-fixed vertices on the boundary as rollers.

    pz : bool / float
        If set to False, no vertical loads apply. If set to True, planar tributary area apply.

    px : bool / float
        If set to False, no x-horizontal loads apply. If set to float will apply a proportion of the vertical pz as px. Example px = 0.2 apply 20% of pz to px.

    px : bool / float
        If set to False, no y-horizontal loads apply. If set to float will apply a proportion of the vertical pz as py. Example py = 0.2 apply 20% of pz to py.

    Returns
    -------
    obj
        FormDiagram.

    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    x_span = x1 - x0
    y_span = y1 - y0
    dx = x_span/division
    dy = y_span/division

    lines = []

    for i in range(division+1):
        for j in range(division+1):
            if i < division and j < division:
                # Vertical Members:
                xa = x0 + dx*i
                ya = y0 + dy*j
                xb = x0 + dx*(i + 1)
                yb = y0 + dy*j
                # Horizontal Members:
                xc = x0 + dx*i
                yc = y0 + dy*j
                xd = x0 + dx*i
                yd = y0 + dy*(j + 1)
                lines.append([[xa,ya,0.0],[xb,yb,0.0]])
                lines.append([[xc,yc,0.0],[xd,yd,0.0]])
                if i == j:
                    # Diagonal Members in + Direction:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*(i + 1)
                    yd = y0 + dy*(j + 1)
                    lines.append([[xc,yc,0.0],[xd,yd,0.0]])
                if i + j == division:
                    # Diagonal Members in - Direction:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*(i - 1)
                    yd = y0 + dy*(j + 1)
                    lines.append([[xc,yc,0.0],[xd,yd,0.0]])
                    if i == (division - 1):
                        xc = x0 + dx*i
                        yc = y0 + dy*j
                        xd = x0 + dx*(i + 1)
                        yd = y0 + dy*(j - 1)
                        lines.append([[xc,yc,0.0],[xd,yd,0.0]])
            else:
                if i == division and j < division:
                    # Vertical Members on last column:
                    xa = x0 + dx*j
                    ya = y0 + dy*i
                    xb = x0 + dx*(j + 1)
                    yb = y0 + dy*i
                    # Horizontal Members:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*i
                    yd = y0 + dy*(j + 1)
                    lines.append([[xa,ya,0.0],[xb,yb,0.0]])
                    lines.append([[xc,yc,0.0],[xd,yd,0.0]])

    form = FormDiagram.from_lines(lines, delete_boundary_face=True)

    form.update_default_vertex_attributes({'is_roller': False})
    form.update_default_vertex_attributes({'is_fixed': False})
    form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form.attributes['loadpath'] = 0
    form.attributes['indset'] = []

    gkey_key = form.gkey_key()

    if fix == 'corners':
        form.set_vertex_attribute(gkey_key[geometric_key([x0,y0,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x0,y1,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x1,y0,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x1,y1,0.0])], 'is_fixed', True)
    else:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            form.set_vertex_attribute(key, 'is_fixed', True)
    
    if rollers:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            if form.get_vertex_attribute(key, 'is_fixed') == False:
                form.set_vertex_attribute(key, 'is_roller', True)

    if pz:
        pzt = 0
        for key in form.vertices():
            form.vertex[key]['pz'] = form.vertex_area(key=key)
            pzt += form.vertex[key]['pz']
        print('Planar load - pzt = {0}'.format(pzt))

    if px:
        pxt = 0
        for key in form.vertices():
            form.vertex[key]['px'] = form.vertex[key]['pz'] * px
            pxt += form.vertex[key]['px']
        print('Load x-direction - pxt = {0}'.format(pxt))

    if py:
        pyt = 0
        for key in form.vertices():
            form.vertex[key]['py'] = form.vertex[key]['pz'] * py
            pyt += form.vertex[key]['py']
        print('Loax y-direction - pyt = {0}'.format(pyt))

    return form

def create_fan_form(xy_span = [[0.0,10.0],[0.0,10.0]], division = 10, fix = 'corners', pz = True, px = None, py = None):

    """ Create a cross form-diagram.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    division: int
        Set the density of the grid in x and y directions.

    fix : string
        Option to select the constrained nodes: 'corners', 'all' are accepted.

    rollers: bool
        Set the non-fixed vertices on the boundary as rollers.

    pz : bool / float
        If set to False, no vertical loads apply. If set to True, planar tributary area apply.

    px : bool / float
        If set to False, no x-horizontal loads apply. If set to float will apply a proportion of the vertical pz as px. Example px = 0.2 apply 20% of pz to px.

    px : bool / float
        If set to False, no y-horizontal loads apply. If set to float will apply a proportion of the vertical pz as py. Example py = 0.2 apply 20% of pz to py.

    Returns
    -------
    obj
        FormDiagram.

    """

    from compas.geometry.transformations.transformations import mirror_point_line
    
    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    x_span = x1 - x0
    y_span = y1 - y0
    xc0 = x0 + x_span/2
    yc0 = y0 + y_span/2
    dx = float(x_span/division)
    dy = float(y_span/division)
    n = int(division/2)
    line_hor = [[x0,yc0,0.0],[xc0,yc0,0.0]]
    line_ver = [[xc0,y0,0.0],[xc0,yc0,0.0]]

    lines = []

    for i in range(n):
        for j in range(n+1):
            # Diagonal Members:
            xa = x0 + dx*i
            ya = y0 + dy*j*i/n
            xb = x0 + dx*(i + 1)
            yb = y0 + dy*j*(i + 1)/n
            lines.append([[xa,ya,0.0],[xb,yb,0.0]])

            a_mirror = mirror_point_line([xa, ya, 0.0], line_hor)
            b_mirror = mirror_point_line([xb, yb, 0.0], line_hor)
            lines.append([a_mirror,b_mirror])
            a_mirror = mirror_point_line(a_mirror, line_ver)
            b_mirror = mirror_point_line(b_mirror, line_ver)
            lines.append([a_mirror,b_mirror])
            a_mirror = mirror_point_line([xa, ya, 0.0], line_ver)
            b_mirror = mirror_point_line([xb, yb, 0.0], line_ver)
            lines.append([a_mirror,b_mirror])

            xa_ = x0 + dx * j * i / n
            ya_ = y0 + dy * i
            xb_ = x0 + dx * j * (i + 1) / n
            yb_ = y0 + dy * (i + 1)
            lines.append([[xa_,ya_,0.0],[xb_,yb_,0.0]])

            a_mirror = mirror_point_line([xa_, ya_, 0.0], line_hor)
            b_mirror = mirror_point_line([xb_, yb_, 0.0], line_hor)
            lines.append([a_mirror,b_mirror])
            a_mirror = mirror_point_line(a_mirror, line_ver)
            b_mirror = mirror_point_line(b_mirror, line_ver)
            lines.append([a_mirror,b_mirror])
            a_mirror = mirror_point_line([xa_, ya_, 0.0], line_ver)
            b_mirror = mirror_point_line([xb_, yb_, 0.0], line_ver)
            lines.append([a_mirror,b_mirror])

            if j < n:
                # Vertical or Horizontal Members:
                xc = x0 + dx * (i + 1)
                yc = y0 + dy * j * (i + 1) / n
                xd = x0 + dx * (i + 1)
                yd = y0 + dy * (j + 1) * (i + 1) / n
                lines.append([[xc,yc,0.0],[xd,yd,0.0]])

                c_mirror = mirror_point_line([xc, yc, 0.0], line_hor)
                d_mirror = mirror_point_line([xd, yd, 0.0], line_hor)
                lines.append([c_mirror, d_mirror])
                c_mirror = mirror_point_line(c_mirror, line_ver)
                d_mirror = mirror_point_line(d_mirror, line_ver)
                lines.append([c_mirror,d_mirror])
                c_mirror = mirror_point_line([xc, yc, 0.0], line_ver)
                d_mirror = mirror_point_line([xd, yd, 0.0], line_ver)
                lines.append([c_mirror, d_mirror])

                xc_ = x0 + dx * j * (i + 1) / n
                yc_ = y0 + dy * (i + 1)
                xd_ = x0 + dx * (j + 1) * (i + 1) / n
                yd_ = y0 + dy * (i + 1)
                lines.append([[xc_,yc_,0.0],[xd_,yd_,0.0]])

                c_mirror = mirror_point_line([xc_, yc_, 0.0], line_hor)
                d_mirror = mirror_point_line([xd_, yd_, 0.0], line_hor)
                lines.append([c_mirror,d_mirror])
                c_mirror = mirror_point_line(c_mirror, line_ver)
                d_mirror = mirror_point_line(d_mirror, line_ver)
                lines.append([c_mirror,d_mirror])
                c_mirror = mirror_point_line([xc_, yc_, 0.0], line_ver)
                d_mirror = mirror_point_line([xd_, yd_, 0.0], line_ver)
                lines.append([c_mirror,d_mirror])


    form = FormDiagram.from_lines(lines)

    form.update_default_vertex_attributes({'is_roller': False})
    form.update_default_vertex_attributes({'is_fixed': False})
    form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form.attributes['loadpath'] = 0
    form.attributes['indset'] = []

    gkey_key = form.gkey_key()

    if fix == 'corners':
        form.set_vertex_attribute(gkey_key[geometric_key([x0,y0,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x0,y1,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x1,y0,0.0])], 'is_fixed', True)
        form.set_vertex_attribute(gkey_key[geometric_key([x1,y1,0.0])], 'is_fixed', True)
    else:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            form.set_vertex_attribute(key, 'is_fixed', True)

    if pz:
        pzt = 0
        for key in form.vertices():
            form.vertex[key]['pz'] = form.vertex_area(key=key)
            pzt += form.vertex[key]['pz']
        print('Planar load - pzt = {0}'.format(pzt))

    if px:
        pxt = 0
        for key in form.vertices():
            form.vertex[key]['px'] = form.vertex[key]['pz'] * px
            pxt += form.vertex[key]['px']
        print('Load x-direction - pxt = {0}'.format(pxt))

    if py:
        pyt = 0
        for key in form.vertices():
            form.vertex[key]['py'] = form.vertex[key]['pz'] * py
            pyt += form.vertex[key]['py']
        print('Loax y-direction - pyt = {0}'.format(pyt))

    return form

def create_dome_form(center = [5.0, 5.0], radius = 5.0, n_radial = 5, n_spikes = 10, r_oculus = 0.0, pz = True, px = None, py = None):

    """ Create a Dome Form-diagram and set common attributes.

    Parameters
    ----------
    center : list
        Planar coordinates of the form-diagram [xc, yc].

    radius: float
        Radius of the form-diagram

    n_radial : int
        Number of meridians on the dome form-diagram.

    n_spikes: int
        Number of spikes from the center.

    r_oculus: float
        Value of the radius of the oculus, if no oculus is present should be set to zero.

    pz : bool / float
        If set to False, no vertical loads apply. If set to True, planar tributary area apply.

    px : bool / float
        If set to False, no x-horizontal loads apply. If set to float will apply a proportion of the vertical pz as px. Example px = 0.2 apply 20% of pz to px.

    px : bool / float
        If set to False, no y-horizontal loads apply. If set to float will apply a proportion of the vertical pz as py. Example py = 0.2 apply 20% of pz to py.

    Returns
    -------
    obj
        FormDiagram.

    """

    xc = center[0]
    yc = center[1]
    theta = 2*math.pi/n_spikes
    r_div = (radius - r_oculus)/n_radial
    lines = []

    from compas_tna.utilities import LoadUpdater

    for nr in range(n_radial+1):
        for nc in range(n_spikes):
            if (r_oculus + nr * r_div) > 0.0:
                # Meridian Elements
                xa = xc + (r_oculus + nr * r_div) * math.cos( theta * nc )
                xb = xc + (r_oculus + nr * r_div) * math.cos( theta * ( nc + 1 ) )
                ya = yc + (r_oculus + nr * r_div) * math.sin( theta * nc)
                yb = yc + (r_oculus + nr * r_div) * math.sin( theta * ( nc + 1 ) )          
                lines.append([[xa, ya, 0.0],[xb, yb, 0.0]])
            
            if nr <= n_radial - 1:
            
                # Radial Elements
                xa = xc + (r_oculus + nr * r_div) * math.cos( theta * nc )
                xb = xc + (r_oculus + (nr + 1) * r_div) * math.cos( theta * nc)
                ya = yc + (r_oculus + nr * r_div) * math.sin( theta * nc)
                yb = yc + (r_oculus + (nr + 1) * r_div) * math.sin( theta * nc)           
                lines.append([[xa, ya, 0.0],[xb, yb, 0.0]])

    form = FormDiagram.from_lines(lines, delete_boundary_edges = True)

    form.update_default_vertex_attributes({'is_roller': False})
    form.update_default_vertex_attributes({'is_fixed': False})
    form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form.attributes['loadpath'] = 0
    form.attributes['indset'] = []

    if r_oculus:
        for key in form.faces():
            centroid = form.face_centroid(key)
            if centroid[0] == xc and centroid[1] == yc:
                form.set_face_attribute(key, 'is_loaded', False)

    [bnds] = form.vertices_on_boundaries()
    for key in bnds:
        form.set_vertex_attribute(key, 'is_fixed', True)

    if pz:
        pzt = 0
        for key in form.vertices():
            form.vertex[key]['pz'] = form.vertex_area(key=key)
            pzt += form.vertex[key]['pz']
        print('Planar load - pzt = {0}'.format(pzt))

    # updateloads = LoadUpdater(form,0)
    # print(updateloads) # HOW UPDATE OPENINGS

    if px:
        pxt = 0
        for key in form.vertices():
            form.vertex[key]['px'] = form.vertex[key]['pz'] * px
            pxt += form.vertex[key]['px']
        print('Load x-direction - pxt = {0}'.format(pxt))

    if py:
        pyt = 0
        for key in form.vertices():
            form.vertex[key]['py'] = form.vertex[key]['pz'] * py
            pyt += form.vertex[key]['py']
        print('Loax y-direction - pyt = {0}'.format(pyt))

    return form

def sym_on_openings(form, xy_span = [[0.0,10.0],[0.0,10.0]], exc_length = 1.0):

    gkey_key = form.gkey_key()

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    bndr = form.vertices_on_boundary()
    lines = [form.edge_coordinates(u,v) for u,v in form.edges()]

    for key in bndr:
        if form.get_vertex_attribute(key, 'is_fixed') == False:
            x, y, _ = form.vertex_coordinates(key)
            if x == x1: 
                lines.append([[x, y, 0.0],[x + exc_length, y, 0.0]])
            if x == x0: 
                lines.append([[x, y, 0.0],[x - exc_length, y, 0.0]])
            if y == y1: 
                lines.append([[x, y, 0.0],[x, y + exc_length, 0.0]])
            if y == y0: 
                lines.append([[x, y, 0.0],[x, y - exc_length, 0.0]])

    form_ = FormDiagram.from_lines(lines, delete_boundary_face = False)
    # form.delete_face(0)

    form_.update_default_vertex_attributes({'is_roller': False})
    form_.update_default_vertex_attributes({'is_fixed': False})
    form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form_.attributes['loadpath'] = 0
    form_.attributes['indset'] = []

    for key_ in form_.vertices():
        coord_ = form_.vertex_coordinates(key_)
        try:
            key = gkey_key[geometric_key(coord_)]
            pz = form.get_vertex_attribute(key, 'pz')
            px = form.get_vertex_attribute(key, 'px')
            py = form.get_vertex_attribute(key, 'py')
            fixed = form.get_vertex_attribute(key, 'is_fixed')
            lb = form.get_vertex_attribute(key, 'lb')
            ub = form.get_vertex_attribute(key, 'ub')
            target = form.get_vertex_attribute(key, 'target')
            form_.set_vertex_attribute(key_, 'pz', value = pz)
            form_.set_vertex_attribute(key_, 'px', value = px)
            form_.set_vertex_attribute(key_, 'py', value = py)
            form_.set_vertex_attribute(key_, 'is_fixed', value = fixed)
            form_.set_vertex_attribute(key_, 'lb', value = lb)
            form_.set_vertex_attribute(key_, 'ub', value = ub)
            form_.set_vertex_attribute(key_, 'target', value = target)
        except:
            form_.set_vertex_attribute(key_, 'is_fixed', value = True)
            form_.set_vertex_attribute(key_, 'pz', value = 0.0)
            form_.set_vertex_attribute(key_, 'px', value = 0.0)
            form_.set_vertex_attribute(key_, 'py', value = 0.0)
            form_.set_vertex_attribute(key_, 'lb', value = 0.0)
            form_.set_vertex_attribute(key_, 'ub', value = 0.0)
            form_.set_vertex_attribute(key_, 'target', value = 0.0)
            ngb = form_.vertex_neighbors(key_)[0]
            form_.set_edge_attribute((key_,ngb), 'is_symmetry', True)

    return form_