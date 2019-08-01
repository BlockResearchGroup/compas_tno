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

from numpy import array
from numpy import float64
from numpy import zeros
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum


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
    'oveview_forces',
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
    lb = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'ub') for key in form.vertices()}
    ub = {geometric_key(form.vertex_coordinates(key)[:2] + [0]): form.get_vertex_attribute(key, 'lb') for key in form.vertices()}
    
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


def oveview_forces(form):

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
            print(qi)
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

    return