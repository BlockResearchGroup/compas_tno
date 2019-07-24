from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from compas.numerical import connectivity_matrix

from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import abs
from numpy import argmin
from numpy import array
from numpy import float64
from numpy import newaxis

from copy import deepcopy

from compas_thrust.algorithms.equilibrium import z_from_form
from compas_thrust.algorithms.equilibrium import update_form
from compas_thrust.algorithms.equilibrium import paralelise_form

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

__all__ =[
    'lagrangian_scale',
    'evaluate_scale',
    'scale_form',
    'scale_fdm',
]

def lagrangian_scale(form):

    # Mapping

    k_i  = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices_where({'is_external': False})]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True, 'is_external' : False})]
    fixed = [k_i[key] for key in form.vertices_where({'is_external': False, 'is_fixed': True})]
    free  = list(set(vertices) - set(fixed))

    m = len(edges)
    n = len(vertices)
    ni = len(free)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    q   = zeros((m, 1))

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            i = k_i[key]
            xyz[i, :] = form.vertex_coordinates(key)
            x[i]  = vertex.get('x')
            y[i]  = vertex.get('y')
            z[i]  = vertex.get('z')
            s[i]  = vertex.get('target')
            px[i] = vertex.get('px', 0)
            py[i] = vertex.get('py', 0)
            pz[i] = vertex.get('pz', 0)

    q   = zeros((m, 1))

    for u,v in form.edges():
        if form.get_edge_attribute((u,v), 'is_edge') == True:
            i = uv_i[(u,v)]
            q[i] = form.get_edge_attribute((u,v), 'q')

    p = pz[free]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Cit = Ci.transpose()
    D = Cit.dot(diags(q.flatten())).dot(C)
    Di = D[:, free]
    Df = D[:, fixed]
    Dt = D.transpose()

    args = q, free, fixed, C, Ci, Cit, Cf, pz, s, z

    A11 = diags([1]*n+[0]).toarray()
    A12 = vstack([Dt.toarray(),-1*p.transpose()])
    A21 = hstack([D.toarray(),-1*p])
    A22 = zeros((ni,ni))
    A = vstack([hstack([A11,A12]),hstack([A21,A22])])
    B = vstack([1*s,zeros((1+ni, 1))])

    sol = spsolve(A, B)
    r = sol[n]
    # r = 0.59
    print('The lagragian scale is: {0:.3f}'.format(r))

    # r = 2.0

    # Update zs

    q = q * r
    Q = diags([q.flatten()], [0])
    D = Cit.dot(Q).dot(C)
    Di = D[:, free]
    Df = D[:, fixed]
    A_= Di
    b_ = pz[free, 0] - Df.dot(z[fixed,0])
    z[free, 0] = spsolve(A_, b_)

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            i = k_i[key]
            [zi] = z[i]
            form.set_vertex_attribute(key, 'z', value = zi)

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True:
            i = uv_i[(u,v)]
            [qi] = q[i]
            form.set_edge_attribute((u,v),'q',value=qi)

    return form

def evaluate_scale(form, function, bounds, n = 100, plot = True):

    r0 = bounds[0]
    stp = (bounds[1]-bounds[0])/n
    x = []
    y = []
    q0 = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]

    form_ = deepcopy(form)

    k_i  = form_.key_index()
    uv_i = form_.uv_index()

    for k in range(n):
        r = r0 + stp * k
        q = q0 * r
        x.append(r)

        for u, v in form_.edges_where({'is_edge': True}):
            i = uv_i[(u,v)]
            [qi] = q[i]
            form_.set_edge_attribute((u,v),'q',value=qi)

        form_ = z_from_form(form_)
        y.append(function(form_))

    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(x, y)

    ax.set(xlabel='Scale (r)', ylabel='Energy (f)',
        title='Evaluate Energy by Scaling')
    ax.grid()

    pos = argmin(y)
    xmin = x[pos]
    ymin = y[pos]
    text= "x={:.3f}, y={:.3f}".format(xmin, ymin)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
            arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmin, ymin), xytext=(0.94,0.96), **kw)

    if plot:
        plt.show()

    return xmin

def scale_fdm(form, r):

    uv_i = form.uv_index()
    q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True, 'is_external' : False})])[:, newaxis]
    q = q * r

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True:
            i = uv_i[(u,v)]
            [qi] = q[i]
            form.set_edge_attribute((u,v),'q',value=qi)

    form = z_from_form(form)

    return form

def scale_form(form,r):

    k_i     = form.key_index()
    uv_i    = form.uv_index()
    vcount  = len(form.vertex)
    anchors = list(form.anchors())
    fixed   = list(form.fixed())
    fixed   = set(anchors + fixed)
    fixed   = [k_i[key] for key in fixed]
    free    = list(set(range(vcount)) - set(fixed))
    edges   = [(k_i[u], k_i[v]) for u, v in form.edges_where({'is_edge': True})]
    xyz     = array(form.get_vertices_attributes('xyz'), dtype=float64)
    p       = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q       = [attr.get('q', 1.0) for u, v, attr in form.edges_where({'is_edge': True}, True)]
    q       = array(q, dtype=float64).reshape((-1, 1))
    C       = connectivity_matrix(edges, 'csr')
    Ci      = C[:, free]
    Cf      = C[:, fixed]
    Cit     = Ci.transpose()
    Ct      = C.transpose()

    q = q * r
    Q = diags([q.ravel()], [0])

    A       = Cit.dot(Q).dot(Ci)
    B       = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A,p[free, 2] - B.dot(xyz[fixed, 2]))

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z']  = xyz[index, 2]

    for u, v, attr in form.edges_where({'is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['q'] = q[index, 0]

    return form
