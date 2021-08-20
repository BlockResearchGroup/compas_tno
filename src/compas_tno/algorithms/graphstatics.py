

from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_force

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax
from compas_tna.utilities import rot90
from compas_tna.utilities import apply_bounds
from compas_tna.utilities import parallelise_sparse

from compas.numerical import connectivity_matrix
from compas.numerical import spsolve_with_known
from compas.numerical import normrow
from compas.numerical import normalizerow
from compas.geometry import distance_point_point_xy
from compas.geometry import angle_vectors_xy

from compas_plotters import MeshPlotter

from numpy import array
from numpy import float64

from scipy.sparse import diags
from scipy.sparse.linalg import factorized

__all__ = [
    'update_tna',
    'force_update_from_form',
    'paralelise_form',
    'initialise_tna',
    'reciprocal_from_form'
]


def update_tna(form, delete_face=True, plots=False, save=False):

    if delete_face:
        form.delete_face(0)

    corners = list(form.vertices_where({'is_fixed': True}))
    print(form)
    form.set_vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
    form.update_boundaries(feet=2)

    # for key in form.edges_where({'_is_external': True}):
    #     form.edge_attribute(key,'q',value=x_reaction)
    #     form.edge_attribute(key,'fmin',value=x_reaction)
    #     form.edge_attribute(key,'fmax',value=x_reaction)
    #     form.edge_attribute(key,'lmin',value=1.00)
    #     form.edge_attribute(key,'lmax',value=1.00)

    for u, v in form.edges_where({'_is_external': False}):
        qi = form.edge_attribute((u, v), 'q')
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a, b)
        form.edge_attribute((u, v), 'fmin', value=qi*lh)
        form.edge_attribute((u, v), 'fmax', value=qi*lh)
        form.edge_attribute((u, v), 'lmin', value=lh)
        form.edge_attribute((u, v), 'lmax', value=lh)

    force = ForceDiagram.from_formdiagram(form)
    horizontal(form, force, display=False)
    # Vertical?

    if plots:
        plot_force(force, form, radius=0.05).show()
        plot_form(form, radius=0.05).show()

    return form, force


def force_update_from_form(force, form):
    """Update the force diagram after modifying the (force densities of) the form diagram.

    Parameters
    ----------
    force : :class:`ForceDiagram`
        The force diagram on which the update is based.
    form : :class:`FormDiagram`
        The form diagram to update.

    Reference
    -------

    See ``compas_tna`` package.

    Returns
    -------
    None
        The form and force diagram are updated in-place.
    """
    # --------------------------------------------------------------------------
    # form diagram
    # --------------------------------------------------------------------------
    vertex_index = form.vertex_index()

    xy = array(form.xy(), dtype=float64)
    edges = [[vertex_index[u], vertex_index[v]] for u, v in form.edges_where({'_is_edge': True})]
    C = connectivity_matrix(edges, 'csr')
    Q = diags([form.q()], [0])
    uv = C.dot(xy)
    # --------------------------------------------------------------------------
    # force diagram
    # --------------------------------------------------------------------------
    _vertex_index = force.vertex_index()

    _known = [_vertex_index[force.anchor()]]
    _xy = array(force.xy(), dtype=float64)
    _edges = force.ordered_edges(form)
    _edges[:] = [(_vertex_index[u], _vertex_index[v]) for u, v in _edges]
    _C = connectivity_matrix(_edges, 'csr')
    _Ct = _C.transpose()
    print(C.shape, _C.shape, Q.shape, uv.shape)
    # --------------------------------------------------------------------------
    # compute reciprocal for given q
    # --------------------------------------------------------------------------
    _xy = spsolve_with_known(_Ct.dot(_C), _Ct.dot(Q).dot(uv), _xy, _known)
    # --------------------------------------------------------------------------
    # update force diagram
    # --------------------------------------------------------------------------
    for vertex, attr in force.vertices(True):
        index = _vertex_index[vertex]
        attr['x'] = _xy[index, 0]
        attr['y'] = _xy[index, 1]


def paralelise_form(form, force, q, alpha=1.0, kmax=100, plot=None, display=False):

    # Update constraints in edges of Form

    uv_i = form.uv_index()

    for u, v in form.edges_where({'_is_edge': True, '_is_external': False}):
        i = uv_i[(u, v)]
        key = (u, v)
        form.edge_attribute(key, 'q', value=q[i])
        # print(q[i])
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a, b)
        f_target = q[i]*lh
        # print(f_target)
        form.edge_attribute(key, 'fmin', value=f_target)
        form.edge_attribute(key, 'fmax', value=f_target)
        # form.edge_attribute(key,'lmin',value=lh)
        # form.edge_attribute(key,'lmax',value=lh)

    if plot:
        plot_form(form).show()
        plot_force(force, form).show()

    # Initialize

    k_i = form.key_index()
    uv_i = form.uv_index()
    vcount = len(form.vertex)
    anchors = list(form.anchors())
    fixed = list(form.fixed())
    fixed = set(anchors + fixed)
    fixed = [k_i[key] for key in fixed]
    free = list(set(range(vcount)) - set(fixed))
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True})]
    xy = array(form.get_vertices_attributes('xy'), dtype=float64)
    lmin = array([attr.get('lmin', 1e-7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    lmax = array([attr.get('lmax', 1e+7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmin = array([attr.get('fmin', 1e-7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmax = array([attr.get('fmax', 1e+7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    C = connectivity_matrix(edges, 'csr')
    Ct = C.transpose()
    CtC = Ct.dot(C)
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cit = Ci.transpose()

    # force = ForceDiagram.from_formdiagram(form)

    _k_i = force.key_index()
    _fixed = list(force.fixed())
    _fixed = [_k_i[key] for key in _fixed]
    _fixed = _fixed or [0]
    _edges = force.ordered_edges(form)
    _xy = array(force.get_vertices_attributes('xy'), dtype=float64)
    _C = connectivity_matrix(_edges, 'csr')
    _Ct = _C.transpose()
    _Ct_C = _Ct.dot(_C)

    _xy[:] = rot90(_xy, +1.0)

    uv = C.dot(xy)
    _uv = _C.dot(_xy)
    length = normrow(uv)
    _length = normrow(_uv)

    t = alpha * normalizerow(uv) + (1 - alpha) * normalizerow(_uv)

    # Paralelize

    for k in range(kmax):
        # apply length bounds
        apply_bounds(length, lmin, lmax)
        apply_bounds(_length, fmin, fmax)
        # print, if allowed
        if display:
            print(k)
        if alpha != 1.0:
            # if emphasis is not entirely on the form
            # update the form diagram
            xy = parallelise_sparse(CtC, Ct.dot(length * t), xy, fixed, 'CtC')
            uv = C.dot(xy)
            length = normrow(uv)
        if alpha != 0.0:
            # if emphasis is not entirely on the force
            # update the force diagram
            _xy = parallelise_sparse(_Ct_C, _Ct.dot(_length * t), _xy, _fixed, '_Ct_C')
            _uv = _C.dot(_xy)
            _length = normrow(_uv)

    f = _length
    q = (f / length).astype(float64)
    q = q.reshape(-1, 1)
    # print('Final qs')
    # print(q)

    _xy[:] = rot90(_xy, -1.0)

    a = [angle_vectors_xy(uv[i], _uv[i], deg=True) for i in range(len(edges))]

    # print(a)

    for key, attr in form.vertices(True):
        i = k_i[key]
        attr['x'] = xy[i, 0]
        attr['y'] = xy[i, 1]
    for u, v, attr in form.edges_where({'_is_edge': True}, True):
        i = uv_i[(u, v)]
        attr['q'] = q[i, 0]
        attr['f'] = f[i, 0]
        attr['_l'] = length[i, 0]
        attr['a'] = a[i]

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    # Update Z

    xyz = array(form.get_vertices_attributes('xyz'), dtype=float64)
    q = q.reshape(-1, 1)
    Q = diags([q.ravel()], [0])
    p = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)

    A = Cit.dot(Q).dot(Ci)
    A_solve = factorized(A)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = A_solve(p[free, 2] - B.dot(xyz[fixed, 2]))

    length = normrow(C.dot(xyz))
    f = q * length
    r = C.transpose().dot(Q).dot(C).dot(xyz) - p

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z'] = xyz[index, 2].tolist()
        attr['_rx'] = r[index, 0]
        attr['_ry'] = r[index, 1]
        attr['_rz'] = r[index, 2]
    for u, v, attr in form.edges_where({'_is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['f'] = f[index, 0]

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    if plot:
        plot_form(form).show()
        plot_force(force, form).show()

    return form, force


def initialise_tna(form, zmax=5.0, method='nodal', plot=False, alpha=100.0, kmax=500, remove_feet=True, display=False):

    corners = list(form.vertices_where({'is_fixed': True}))
    form.vertices_attribute('is_anchor', True, keys=corners)
    form.edges_attribute('fmin', 0.0)
    form.edges_attribute('fmax', 10.0)
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), '_is_edge') is False:
            leaves = True
            break
    if leaves is False:
        form.update_boundaries()

    force = ForceDiagram.from_formdiagram(form)
    if plot:
        print('Plot of Primal')
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
        plotter.draw_vertices(radius=0.05)
        plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_anchor': True})], radius=0.10, facecolor='000000')
        plotter.show()
        print('Plot of Dual')
        force.plot()

    if method == 'nodal':
        horizontal_nodal(form, force, alpha=alpha, kmax=kmax)
    else:
        horizontal(form, force, alpha=alpha, kmax=kmax)

    vertical_from_zmax(form, zmax)
    # if leaves is False and remove_feet is True:
    #     form = form.remove_feet()

    if plot:
        print('Plot of Reciprocal')
        force.plot()

    return form


def reciprocal_from_form(form, zmax=5.0, method='algebraic', plot=False, alpha=100.0, kmax=500, remove_feet=True, display=False):

    corners = list(form.vertices_where({'is_fixed': True}))
    form.vertices_attribute('is_anchor', True, keys=corners)

    for key in form.vertices_where({'is_fixed': True}):
        break
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), '_is_edge') is False:
            leaves = True
            break
    if leaves is False:
        form.update_boundaries()

    for u, v in form.edges():
        if form.edge_attribute((u, v), '_is_external') is False:
            qi = form.edge_attribute((u, v), 'q')
            a = form.vertex_coordinates(u)
            b = form.vertex_coordinates(v)
            lh = distance_point_point_xy(a, b)
            form.edge_attribute((u, v), 'fmin', value=qi*lh)
            form.edge_attribute((u, v), 'fmax', value=qi*lh)
            form.edge_attribute((u, v), 'lmin', value=lh)
            form.edge_attribute((u, v), 'lmax', value=lh)

    if plot:
        plot_form(form).show()

    force = ForceDiagram.from_formdiagram(form)
    if plot:
        print('Plot of Primal')
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_vertices(radius=0.05)
        plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.10, facecolor='000000')
        plotter.show()
        print('Plot of Dual')
        force.plot()

    if method == 'algebraic':
        force_update_from_form(force, form)
    else:
        horizontal_nodal(form, force, alpha=alpha, kmax=kmax)

    # Vertical Equilibrium with no updated loads

    if leaves is False and remove_feet is True:
        form = form.remove_feet()

    if plot:
        print('Plot of Reciprocal')
        force.plot()

    return form, force
