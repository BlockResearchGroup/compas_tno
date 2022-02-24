from compas_tno.diagrams import ForceDiagram

from compas_tna.equilibrium import horizontal_nodal

from compas_tno.algorithms.equilibrium import vertical_equilibrium_fdm
from compas_tno.plotters import TNOPlotter

from compas.numerical import connectivity_matrix
from compas.numerical import spsolve_with_known

from numpy import array
from numpy import float64

from scipy.sparse import diags


def form_update_with_parallelisation(form, zmax=5.0, plot=False, alpha=100.0, kmax=500):
    """Parallelise a TNO form diagram using as starting point the centroidal dual as ForceDiagram.

    Parameters
    ----------
    form : :class:`FormDiagram`
        The form diagram to update.
    zmax : float, optional
        The maximum height of a node in the thrust network.
        The default value is ``5.0``.
    plot : bool, optional
        If plots with form and force should display on the screen.
        The default value is ``False``.
    alpha : float, optional
        Weighting factor for computation of the target vectors (the default is
        100.0, which implies that the target vectors are the edges of the form diagram).
        If 0.0, the target vectors are the edges of the force diagram.
    kmax : int, optional
       Maximum number of iterations.
       Default is ``500``.

    Reference
    ---------
    See ``compas_tna`` package.

    Returns
    -------
    None
        The form and force diagram are updated in-place.
    """
    # mark supports as 'is_anchor'
    corners = list(form.vertices_where({'is_fixed': True}))
    form.vertices_attribute('is_anchor', True, keys=corners)

    # update bounds adding a new face to the open edges along boundaries
    # does not apply for continuous boundary cases, for such, the edges
    # on the boundary must be previously set with '_is_edge': False
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), '_is_edge') is False:
            leaves = True
            break
    if leaves is False:
        form.update_boundaries()
    force = ForceDiagram.from_formdiagram(form)

    # plot of diagrams @ initial state
    if plot:
        plotter = TNOPlotter(form=form, force=force, figsize=(16, 8))
        plotter.draw_mesh()
        plotter.draw_force()
        plotter.show()

    # horizontal equilibrium
    horizontal_nodal(form, force, alpha=alpha, kmax=kmax)

    # ****** change compression to negative (TNO convention)
    q = [form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]
    for index, edge in enumerate(form.edges_where({'_is_edge': True})):
        form.edge_attribute(edge, 'q', -1 * q[index])

    # update vertical equilibrium (only z coordinates change)
    vertical_equilibrium_fdm(form, zmax=zmax)

    # plot of diagrams @ final state
    if plot:
        plotter = TNOPlotter(form=form, force=force, figsize=(16, 8))
        plotter.draw_mesh()
        plotter.draw_force()
        plotter.show()

    return


def reciprocal_from_form(form, plot=False):
    """Find a reciprocal form diagram from the stored force densities in thisdiagram
    This is an algebraic solution to the duality. It should work as long as the force densities
    in the form diagram are equilibrated.

    Parameters
    ----------
    form : :class:`FormDiagram`
        The form diagram to update.
    plot : bool, optional
        If plots with form and force should display on the screen.
        The default value is ``False``.

    Reference
    ---------
    See ``compas_ags`` package.

    Returns
    -------
    force : :class:`ForceDiagram`
        The reciprocal force diagram.

    """
    # mark supports as 'is_anchor'
    corners = list(form.vertices_where({'is_fixed': True}))
    form.vertices_attribute('is_anchor', True, keys=corners)

    # update bounds adding a new face to the open edges along boundaries
    # does not apply for continuous boundary cases, for such, the edges
    # on the boundary must be previously set with '_is_edge': False
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), '_is_edge') is False:
            leaves = True
            break
    if leaves is False:
        init_faces = list(form.faces())
        form.update_boundaries()


    force = ForceDiagram.from_formdiagram(form)

    if plot:
        plotter = TNOPlotter(form=form, force=force, figsize=(16, 8))
        plotter.draw_mesh()
        plotter.draw_force()
        plotter.show()

    force_update_from_form(force, form)

    if leaves is False:
        end_faces = list(form.faces())
        new_faces = list(set(end_faces) - set(init_faces))
        for fkey in new_faces:
            form.delete_face(fkey)

    if plot:
        print('Plot of Reciprocal')
        force.plot()

    return force


def force_update_from_form(force, form):
    """Update the force diagram after modifying the (force densities of) the form diagram.
    This is an algebraic solution to the duality. It works as long as the force densities
    in the form diagram are equilibrated.

    Parameters
    ----------
    force : :class:`ForceDiagram`
        The force diagram on which the update is based.
    form : :class:`FormDiagram`
        The form diagram to update.

    Reference
    ---------

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
    # --------------------------------------------------------------------------
    # compute reciprocal for given q
    # --------------------------------------------------------------------------
    _xy = spsolve_with_known(_Ct.dot(_C), _Ct.dot(Q).dot(uv), _xy, _known)
    # --------------------------------------------------------------------------
    # rotate force diagram to make it parallel to the form diagram
    # use CCW direction (opposite of cycle direction)
    # --------------------------------------------------------------------------
    _x, _y = zip(*_xy)
    _xy[:] = [list(item) for item in zip([-_ for _ in _y], _x)]
    # --------------------------------------------------------------------------
    # update force diagram
    # --------------------------------------------------------------------------
    for vertex, attr in force.vertices(True):
        index = _vertex_index[vertex]
        attr['x'] = _xy[index, 0]
        attr['y'] = _xy[index, 1]
        attr['z'] = 0.0
