from compas_tno.diagrams import ForceDiagram

from compas_tna.equilibrium import horizontal_nodal

from compas_tno.algorithms.equilibrium import vertical_equilibrium_fdm
from compas_tno.plotters import TNOPlotter

from compas.numerical import connectivity_matrix
from compas.numerical import spsolve_with_known

from numpy import array
from numpy import float64

from scipy.sparse import diags


def form_update_with_parallelisation(form, zmax=5.0, plot=False, printout=False, alpha=100.0, kmax=500, callback=None):
    """Parallelise a TNO form diagram using as starting point the centroidal dual as ForceDiagram.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.
    zmax : float, optional
        The maximum height of a node in the thrust network.
        The default value is ``5.0``.
    plot : bool, optional
        If plots with form and force should display on the screen.
        The default value is ``False``.
    plot : bool, optional
        If prints should appear in the screen with key information.
        The default value is ``False``.
    alpha : float, optional
        Weighting factor for computation of the target vectors (the default is
        100.0, which implies that the target vectors are the edges of the form diagram).
        If 0.0, the target vectors are the edges of the force diagram.
    kmax : int, optional
       Maximum number of iterations.
       Default is ``500``.
    callback : callable, optional
        A callback function to be called at every iteration of the parallelisation algorithm.
        The callback should take the current iterand, the coordinates of the form diagram,
        and the coordinates of the force diagram as input parameters.
        Default is ``None``.

    Nores
    ---------
    See `compas tna <https://blockresearchgroup.github.io/compas_tna/>`_  package.

    Returns
    -------
    force: ForceDiagram
        The force diagram after the process.
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

    # ad = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/force-0.json'
    # force.to_json(ad)

    # for force_edge in force.edges():
    #     form_edge = force.dual_edge(force_edge)
    #     q = form.edge_attribute(form_edge, 'q')
    #     p1, p2 = form.edge_coordinates(*form_edge)
    #     lh = distance_point_point_xy(p1, p2)
    #     lmin = lmax = abs(q * lh)/10
    #     force.edge_attribute(force_edge, 'lmax', lmax)
    #     # if form.edge_attribute(form_edge, 'is_ind'):
    #     #     force.edge_attribute(force_edge, 'lmin', lmin)

    # def callback(k, xyz_force, edges):
    #     print('Iteration:', k)
    #     for i, key in enumerate(force.vertices()):
    #         force.vertex_attribute(key, 'x', xyz_force[i][0])
    #         force.vertex_attribute(key, 'y', xyz_force[i][1])

    #     # ad = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/force-' + str(k + 1) + '.json'
    #     # force.to_json(ad)

    #     plotter = TNOPlotter(force=force)
    #     plotter.draw_force()
    #     plotter.show()

    # plot of diagrams @ initial state
    if plot:
        plotter = TNOPlotter(form=form, force=force, figsize=(16, 8))
        plotter.draw_mesh()
        plotter.draw_force()
        plotter.show()

    # horizontal equilibrium
    horizontal_nodal(form, force, alpha=alpha, kmax=kmax, callback=callback)

    if printout:
        dev = form.edges_attribute('_a')
        print('max/min deviations:', max(dev), min(dev))

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

    return force


def reciprocal_from_form(form, plot=False, restore_form_topology=True):
    """Find a reciprocal form diagram from the stored force densities in thisdiagram
    This is an algebraic solution to the duality. It should work as long as the force densities
    in the form diagram are equilibrated.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.
    plot : bool, optional
        If plots with form and force should display on the screen.
        The default value is ``False``.
    restore_form_topology : bool, optional
        If form topology should be restored.
        The default value is ``True``.

    Notes
    ---------
    See `compas ags <https://blockresearchgroup.github.io/compas_ags/>`_  package.

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

    # print('edges', form.number_of_edges())
    # print('vertices', form.number_of_vertices())
    # print('faces', form.number_of_faces())

    force = ForceDiagram.from_formdiagram(form)

    # print('force edges', force.number_of_edges())
    # print('force vertices', force.number_of_vertices())
    # print('force faces', force.number_of_faces())

    if plot:
        plotter = TNOPlotter(form=form, force=force, figsize=(16, 8))
        plotter.draw_mesh()
        plotter.draw_force()
        plotter.show()

    force_update_from_form(force, form)

    # this modifies "back" the form diagram in case open edges in the boundary are present
    # It will delete the added faces above.
    if restore_form_topology:
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
    force : :class:`~compas_tno.diagrams.ForceDiagram`
        The force diagram on which the update is based.
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.

    Notes
    ---------
    See `compas tna <https://blockresearchgroup.github.io/compas_tna/>`_ package.

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
    # rotate the force diagram 90 degrees in CW direction
    # this way the relation between the two diagrams is easier to read
    # --------------------------------------------------------------------------
    # _x, _y = zip(*_xy)
    # _xy[:] = [list(item) for item in zip(_y, [-_ for _ in _x])]
    # --------------------------------------------------------------------------
    # update force diagram
    # --------------------------------------------------------------------------
    for vertex, attr in force.vertices(True):
        index = _vertex_index[vertex]
        attr['x'] = _xy[index, 0]
        attr['y'] = _xy[index, 1]
        attr['z'] = 0.0
