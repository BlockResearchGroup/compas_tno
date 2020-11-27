from numpy import zeros
import math
from compas_tno.datastructures import MeshDos
from compas.geometry import norm_vector

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

        if intrados.vertex_attribute(key, '_is_outside'):
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
        if intrados.vertex_attribute(key, '_is_outside'):
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
