from numpy import zeros

from compas.geometry import norm_vector

import math


__all__ = [
    'general_ub_lb_update_with_t_middle_constant',
    'general_db_with_t_middle_constant',
    'general_ub_lb_update_with_t_middle_variable',
    'general_db_with_t_middle_variable',
    'general_ub_lb_update_with_t_intrados',
    'general_db_with_t_intrados',
    'general_ub_lb_update_with_s',
    'general_dub_dlb_with_s',
    'general_ub_lb_update_with_n',
    'general_dub_dlb_with_n'
    ]


def general_ub_lb_update_with_t_middle_constant(thk, s, middle, t):

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
    """ This general  update takes into consideration nub and nlb. scaled, so t is a (%)
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

    dub = zeros((len(lb), 1))
    dlb = zeros((len(lb), 1))

    i = 0
    for key in intrados.vertices():
        normal = intrados.vertex_attribute(key, 'n')
        dev = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))  # 1/cos(a)
        dub[i] = 1 * dev * norm_vector(normal)  # Experimenting with this normal norm!
        i += 1

    return dub, dlb


def general_ub_lb_update_with_s(ub, lb, s):  # s represents the "half-portion" of the section that is still remaining

    ub_update = zeros((len(ub), 1))
    lb_update = zeros((len(ub), 1))

    for i in range(len(ub)):
        ub_update[i] = ub[i] - (ub[i] - lb[i]) * s
        lb_update[i] = lb[i] + (ub[i] - lb[i]) * s

    return ub_update, lb_update


def general_dub_dlb_with_s(ub, lb):

    dub = zeros((len(ub), 1))
    dlb = zeros((len(ub), 1))

    for i in range(len(ub)):
        dub[i] = - (ub[i] - lb[i])
        dlb[i] = + (ub[i] - lb[i])

    return dub, dlb


def general_ub_lb_update_with_n(ub, lb, n, intrados, extrados, t):  # n represents the magnitude of the normal offset vector

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
