from compas_tno.shapes.dome import dome_zt_update
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update
from compas_tno.utilities.interpolation import get_shape_middle_pattern

from numpy import zeros


__all__ = [
    'compute_form_initial_lengths',
    'compute_edge_stiffness'
]


def compute_form_initial_lengths(form):
    """Compute the length of each edge based on the projection of the pattern onto the middle surface.
    Note: the middle surface is obtained based on the "target" attributes of the form diagram that should
    be set before invoking this method.

    Parameters
    ----------
    form : ::FormDiagram::
        Form diagram to compute the lengths.

    Returns
    -------
    lengths : array (m x 1)
        A vector (m x 1) containing the calculated edge lengths.
    """

    form_ = form.copy()

    m = form.number_of_real_edges()
    lengths = zeros((m, 1))

    for key in form_.vertices():
        z = form.vertex_attribute(key, 'target')
        form_.vertex_attribute(key, 'z', value=z)

    i = 0
    for edge in form.edges():
        lengths[i] = form_.edge_length(*edge)
        form.edge_attribute(edge, 'l0', lengths[i])
        i += 1

    return lengths


def compute_edge_stiffness(form, lengths=None, area=0.125, E=20E+6):
    """Compute the stiffness of each edge based on the initial lengths, constant area and Young Modulus E.

    Parameters
    ----------
    form : ::FormDiagram::
        Form diagram to compute the stiffness.
    lengths : array (m x 1)
        The lengths considered in each edge.
        The default value is ``None`` in which case the initial lengths are calculated based on the middle surface.
    area : float
        The constant area applied to each edge.
        The default value is ``0.125`` which corresponds to a square section of 0.25 x 0.25 m.
    E : float
        The young modulus of the material.
        The default value is ``20E+6`` which corresponds to a 20 GPa commonly used for concrete.

    Returns
    -------
    k : array (m x 1)
        The stiffness of each edge.

    Note
    -------
    The stiffness is based on a the fo .
    """

    if lengths is None:
        lengths = compute_form_initial_lengths(form)
    m = len(lengths)
    k = zeros((m, 1))

    for index, edge in enumerate(list(form.edges_where({'_is_edge': True}))):
        k[index] = E * area / lengths[index]
        form.edge_attribute(edge, 'k', k[index])

    return k
