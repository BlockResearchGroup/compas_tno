from numpy import array
from numpy import float64
from scipy.sparse import diags

from compas.linalg import spsolve_with_known
from compas.matrices import connectivity_matrix
from compas_tna.equilibrium import horizontal_nodal
from compas_tno.algorithms.equilibrium import vertical_equilibrium_fdm
from compas_tno.diagrams import ForceDiagram
from compas_tno.diagrams import FormDiagram


def form_update_with_parallelisation(
    form: FormDiagram,
    zmax: float = 5.0,
    printout: bool = False,
    alpha: float = 100.0,
    kmax: int = 500,
    callback: callable = None,
) -> ForceDiagram:
    """Parallelise a TNO form diagram using as starting point the centroidal dual as ForceDiagram.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.
    zmax : float, optional
        The maximum height of a node in the thrust network.
        The default value is ``5.0``.
    printout : bool, optional
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

    Returns
    -------
    :class:`ForceDiagram`
        The force diagram after the process.

    Notes
    -----
    See `compas tna <https://blockresearchgroup.github.io/compas_tna/>`_  package.

    """
    # mark supports as 'is_anchor'
    corners = list(form.vertices_where({"": True}))
    form.vertices_attribute("is_anchor", True, keys=corners)

    # update bounds adding a new face to the open edges along boundaries
    # does not apply for continuous boundary cases, for such, the edges
    # on the boundary must be previously set with '_is_edge': False
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), "_is_edge") is False:
            leaves = True
            break
    if leaves is False:
        form.update_boundaries()
    force = ForceDiagram.from_formdiagram(form)

    # horizontal equilibrium
    horizontal_nodal(form, force, alpha=alpha, kmax=kmax, callback=callback)

    if printout:
        dev = form.edges_attribute("_a")
        print("max/min deviations:", max(dev), min(dev))

    # change compression to negative (TNO convention)
    q = [form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]
    for index, edge in enumerate(form.edges_where({"_is_edge": True})):
        form.edge_attribute(edge, "q", -1 * q[index])

    # update vertical equilibrium (only z coordinates change)
    vertical_equilibrium_fdm(form, zmax=zmax)

    return force


def reciprocal_from_form(
    form: FormDiagram,
    restore_form_topology: bool = True,
) -> ForceDiagram:
    """Find a reciprocal form diagram from the stored force densities in thisdiagram
    This is an algebraic solution to the duality. It should work as long as the force densities
    in the form diagram are equilibrated.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.
    restore_form_topology : bool, optional
        If form topology should be restored.
        The default value is ``True``.

    Returns
    -------
    :class:`ForceDiagram`
        The reciprocal force diagram.

    Notes
    -----
    See `compas ags <https://blockresearchgroup.github.io/compas_ags/>`_  package.

    """
    # mark supports as 'is_anchor'
    corners = list(form.vertices_where({"": True}))
    form.vertices_attribute("is_anchor", True, keys=corners)

    # update bounds adding a new face to the open edges along boundaries
    # does not apply for continuous boundary cases, for such, the edges
    # on the boundary must be previously set with '_is_edge': False
    leaves = False
    for u, v in form.edges_on_boundary():
        if form.edge_attribute((u, v), "_is_edge") is False:
            leaves = True
            break
    if leaves is False:
        init_faces = list(form.faces())
        form.update_boundaries()

    force = ForceDiagram.from_formdiagram(form)

    force_update_from_form(force, form)

    # this modifies "back" the form diagram in case open edges in the boundary are present
    # It will delete the added faces above.
    if restore_form_topology:
        if leaves is False:
            end_faces = list(form.faces())
            new_faces = list(set(end_faces) - set(init_faces))
            for fkey in new_faces:
                form.delete_face(fkey)

    return force


def force_update_from_form(force: ForceDiagram, form: FormDiagram) -> None:
    """Update the force diagram after modifying the (force densities of) the form diagram.
    This is an algebraic solution to the duality. It works as long as the force densities
    in the form diagram are equilibrated.

    Parameters
    ----------
    force : :class:`~compas_tno.diagrams.ForceDiagram`
        The force diagram on which the update is based.
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to update.

    Returns
    -------
    None
        The form and force diagram are updated in-place.

    Notes
    -----
    See `compas tna <https://blockresearchgroup.github.io/compas_tna/>`_ package.

    """
    # --------------------------------------------------------------------------
    # form diagram
    # --------------------------------------------------------------------------
    vertex_index = form.vertex_index()

    xy = array(form.xy(), dtype=float64)
    edges = [[vertex_index[u], vertex_index[v]] for u, v in form.edges_where({"_is_edge": True})]
    C = connectivity_matrix(edges, "csr")
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
    _C = connectivity_matrix(_edges, "csr")
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
        attr["x"] = _xy[index, 0]
        attr["y"] = _xy[index, 1]
        attr["z"] = 0.0
