from numpy import zeros

from compas.geometry import norm_vector

import math


def general_ub_lb_update_with_t_middle_constant(thk, s, middle, t):
    """Update upper and lower bounds for a general surface for a given middle surface.

    Parameters
    ----------
    thk : float
        Thickness of the structure
    s : float
        Middle surface heights
    middle : :class:`~compas_tno.shapes.MeshDos`
        Middle surface mesh
    t : float
        Value to assumed to vertices without intrados projection

    Returns
    -------
    ub_update : array
        Upper bound limits updated
    lb_update : array
        Lower bound limits updated
    """

    ub_update = zeros((len(s), 1))
    lb_update = zeros((len(s), 1))

    i = 0
    for key in middle.vertices():
        normal = middle.vertex_attribute(key, 'n')
        dev = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))  # 1/cos(a)
        ub_update[i] = s[i] + thk/2 * dev * norm_vector(normal)  # Experimenting with this normal norm!
        if middle.vertex_attribute(key, 'is_outside'):
            lb_update[i] = t
        else:
            lb_update[i] = s[i] - thk/2 * dev * norm_vector(normal)
        i += 1

    return ub_update, lb_update


def general_db_with_t_middle_constant(s, middle):
    """Update sensitivity in the upper and lower bounds for a general surface for a given middle surface.

    Parameters
    ----------
    s : float
        Middle surface heights
    middle : MeshDos
        Middle surface mesh

    Returns
    -------
    dub : array
        Sensitivities of the upper bound
    dlb : array
        Sensitivities of the lower bound
    """

    dub = zeros((len(s), 1))
    dlb = zeros((len(s), 1))

    i = 0
    for key in middle.vertices():
        normal = middle.vertex_attribute(key, 'n')
        dev = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))  # 1/cos(a)
        dub[i] = 1/2 * dev * norm_vector(normal)  # Experimenting with this normal norm!
        if middle.vertex_attribute(key, 'is_outside'):
            dlb[i] = 0.0
        else:
            dlb[i] = - 1/2 * dev * norm_vector(normal)
        i += 1

    return dub, dlb


def general_ub_lb_update_with_t_middle_variable(thk_alfa, s, middle, t):
    """Update upper and lower bounds for a general surface for a given middle surface. This general  update takes into consideration nub and nlb. scaled, so t is a (%)

    Parameters
    ----------
    thk_alfa : float
        Thickness of the structure as a percentage of nub, nlb
    s : float
        Middle surface heights
    middle : :class:`~compas_tno.shapes.MeshDos`
        Middle surface mesh
    t : float
        Value to assumed to vertices without intrados projection

    Returns
    -------
    ub_update : array
        Upper bound limits updated
    lb_update : array
        Lower bound limits updated
    """

    ub_update = zeros((len(s), 1))
    lb_update = zeros((len(s), 1))

    i = 0
    for key in middle.vertices():
        # normal = middle.vertex_attribute(key, 'n')
        nub = middle.vertex_attribute(key, 'nub')
        nlb = middle.vertex_attribute(key, 'nlb')
        dev_ub = 1/math.sqrt(1/(1 + (nub[0]**2 + nub[1]**2)/nub[2]**2))  # 1/cos(a)
        if middle.vertex_attribute(key, 'is_outside'):
            dev_lb = 0.0
        else:
            dev_lb = 1/math.sqrt(1/(1 + (nlb[0]**2 + nlb[1]**2)/nlb[2]**2))  # 1/cos(a)
        ub_update[i] = + thk_alfa * dev_ub * norm_vector(nub)  # Experimenting with this normal norm!
        lb_update[i] = - thk_alfa * dev_lb * norm_vector(nlb)  # check if neeed to divide by 2
        # if intrados.vertex_attribute(key, 'is_outside'):
        #     lb_update[i] = t
        i += 1

    return ub_update + s, lb_update + s


def general_db_with_t_middle_variable(s, middle):
    """Update sensitivity in the upper and lower bounds for a general surface for a given middle surface considering the variable as a percentage.

    Parameters
    ----------
    s : float
        Middle surface heights
    middle : :class:`~compas_tno.shapes.MeshDos`
        Middle surface mesh

    Returns
    -------
    dub : array
        Sensitivities of the upper bound
    dlb : array
        Sensitivities of the lower bound
    """

    dub = zeros((len(s), 1))
    dlb = zeros((len(s), 1))

    i = 0
    for key in middle.vertices():
        # normal = middle.vertex_attribute(key, 'n')
        nub = middle.vertex_attribute(key, 'nub')
        nlb = middle.vertex_attribute(key, 'nlb')
        dev_ub = 1/math.sqrt(1/(1 + (nub[0]**2 + nub[1]**2)/nub[2]**2))  # 1/cos(a)
        if middle.vertex_attribute(key, 'is_outside'):
            dev_lb = 0.0
        else:
            dev_lb = 1/math.sqrt(1/(1 + (nlb[0]**2 + nlb[1]**2)/nlb[2]**2))  # 1/cos(a)
        dub[i] = + 1 * dev_ub * norm_vector(nub)
        dlb[i] = - 1 * dev_lb * norm_vector(nlb)
        i += 1

    return dub, dlb


def general_ub_lb_update_with_t_intrados(thk, lb, intrados, t):
    """Update upper and lower bounds for a general surface for a given lower bound (intrados) surface.

    Parameters
    ----------
    thk_alfa : float
        Thickness of the structure as a percentage of nub, nlb
    s : float
        Middle surface heights
    middle : :class:`~compas_tno.shapes.MeshDos`
        Middle surface mesh
    t : float
        Value to assumed to vertices without intrados projection

    Returns
    -------
    ub_update : array
        Upper bound limits updated
    lb_update : array
        Lower bound limits updated
    """

    ub_update = zeros((len(lb), 1))
    lb_update = lb

    i = 0
    for key in intrados.vertices():
        normal = intrados.vertex_attribute(key, 'n')
        dev = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))  # 1/cos(a)
        ub_update[i] = lb[i] + 1.0 * thk * dev * norm_vector(normal)  # Experimenting with this normal norm!
        i += 1

    return ub_update, lb_update


def general_db_with_t_intrados(lb, intrados):
    """Update sensitivity in the upper and lower bounds for a general surface for a given lower bound (intrados) surface.

    Parameters
    ----------
    s : float
        Middle surface heights
    middle : :class:`~compas_tno.shapes.MeshDos`
        Middle surface mesh

    Returns
    -------
    dub : array
        Sensitivities of the upper bound
    dlb : array
        Sensitivities of the lower bound
    """

    dub = zeros((len(lb), 1))
    dlb = zeros((len(lb), 1))

    i = 0
    for key in intrados.vertices():
        normal = intrados.vertex_attribute(key, 'n')
        dev = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))  # 1/cos(a)
        dub[i] = 1 * dev * norm_vector(normal)  # Experimenting with this normal norm!
        i += 1

    return dub, dlb


def general_ub_lb_update_with_s(ub, lb, s):
    """Update upper and lower bounds for a general surface based on upper and lower bounds that are squized by the ``s`` variable.

    Parameters
    ----------
    ub : array
        Upper bound original limits
    lb : array
        Lower bound original limits
    s : float
        Parameter to reduce the thickness (as a percentual of the original section)

    Returns
    -------
    ub_update : array
        Upper bound limits updated
    lb_update : array
        Lower bound limits updated
    """

    ub_update = zeros((len(ub), 1))
    lb_update = zeros((len(ub), 1))

    for i in range(len(ub)):
        ub_update[i] = ub[i] - (ub[i] - lb[i]) * s
        lb_update[i] = lb[i] + (ub[i] - lb[i]) * s

    return ub_update, lb_update


def general_dub_dlb_with_s(ub, lb):
    """Update sensitivity in the upper and lower bounds for a general surface for a given lower bound (intrados) surface.

    Parameters
    ----------
    ub : array
        Current bound limits
    lb : array
        Current bound limits

    Returns
    -------
    dub : array
        Sensitivities of the upper bound
    dlb : array
        Sensitivities of the lower bound
    """

    dub = zeros((len(ub), 1))
    dlb = zeros((len(ub), 1))

    for i in range(len(ub)):
        dub[i] = - (ub[i] - lb[i])
        dlb[i] = + (ub[i] - lb[i])

    return dub, dlb


def general_ub_lb_update_with_n(ub, lb, n, intrados, extrados, t):
    """Update upper and lower bounds for a general surface for two given upper and lower bound surfaces.

    Parameters
    ----------
    ub : array
        Current bound limits
    lb : array
        Current bound limits
    n : float
        Magnitude of the offset from the original bounds
    intrados : MeshDos
        Intrados surface mesh
    extrados : MeshDos
        Extrados surface mesh
    t : float
        Value to assumed to vertices without intrados projection

    Returns
    -------
    ub_update : array
        Upper bound limits updated
    lb_update : array
        Lower bound limits updated
    """

    ub_update = zeros((len(ub), 1))
    lb_update = zeros((len(ub), 1))

    # Intrados and Extrados must me diagrams with exactly same topology as the form diagram: i.e. projections of the form diagram with z = ub,lb.

    i = 0
    for key in intrados.vertices():
        normal_intra = intrados.vertex_attribute(key, 'n')
        normal_extra = extrados.vertex_attribute(key, 'n')

        deviation_extra = 1/math.sqrt(1/(1 + (normal_extra[0]**2 + normal_extra[1]**2)/normal_extra[2]**2))  # 1/cos(a)
        ub_update[i] = ub[i] - n * deviation_extra * norm_vector(normal_extra)  # Experimenting with this normal norm!

        if intrados.vertex_attribute(key, 'is_outside'):
            lb_update[i] = t
        else:
            deviation_intra = 1/math.sqrt(1/(1 + (normal_intra[0]**2 + normal_intra[1]**2)/normal_intra[2]**2))  # 1/cos(a)
            lb_update[i] = lb[i] + n * deviation_intra * norm_vector(normal_intra)  # Experimenting with this normal norm!
        i += 1

    return ub_update, lb_update


def general_dub_dlb_with_n(ub, lb, n, intrados, extrados, t):
    """Sensitivity upper and lower bounds for a general surface for two given upper and lower bound surfaces.

    Parameters
    ----------
    ub : array
        Current bound limits
    lb : array
        Current bound limits
    n : float
        Magnitude of the offset from the original bounds
    intrados : :class:`~compas_tno.shapes.MeshDos`
        Intrados surface mesh
    extrados : :class:`~compas_tno.shapes.MeshDos`
        Extrados surface mesh
    t : float
        Value to assumed to vertices without intrados projection

    Returns
    -------
    dub : array
        Sensitivity of upper bound limits
    dlb : array
        Sensitivity of lower bound limits
    """

    dub = zeros((len(ub), 1))
    dlb = zeros((len(ub), 1))

    i = 0
    for key in intrados.vertices():
        normal_intra = intrados.vertex_attribute(key, 'n')
        normal_extra = extrados.vertex_attribute(key, 'n')
        x, y, z_lb = intrados.vertex_coordinates(key)
        x, y, z_ub = extrados.vertex_coordinates(key)
        if intrados.vertex_attribute(key, 'is_outside'):
            deviation_intra = 0.0
        else:
            deviation_intra = 1/math.sqrt(1/(1 + (normal_intra[0]**2 + normal_intra[1]**2)/normal_intra[2]**2))  # 1/cos(a)
        deviation_extra = 1/math.sqrt(1/(1 + (normal_extra[0]**2 + normal_extra[1]**2)/normal_extra[2]**2))  # 1/cos(a)
        dlb[i] = + deviation_intra * norm_vector(normal_intra)  # Experimenting with this normal norm!
        dub[i] = - deviation_extra * norm_vector(normal_extra)  # Experimenting with this normal norm!
        i += 1

    return dub, dlb


def _invert_vector(vector):
    return [-1*vector[0], -1*vector[1], -1*vector[2]]


def general_b_update_with_n(b, n, fixed):
    return


def general_db_with_n(b, n, fixed):

    return
