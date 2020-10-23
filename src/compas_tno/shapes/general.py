from numpy import zeros
import math
from compas_tno.datastructures import MeshDos

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


def general_ub_lb_update_with_n(ub, lb, n, intrados, extrados, t):  # n represents the magnitude of the normal vector

    ub_update = zeros((len(ub), 1))
    lb_update = zeros((len(ub), 1))

    # Intrados and Extrados must me diagrams with exactly same topology as the form diagram: i.e. projections of the form diagram with z = ub,lb.

    # vertices, faces = intrados.to_vertices_and_faces()
    # intrados = intrados.copy()
    # extrados = extrados.copy()
    # vertices, faces = extrados.to_vertices_and_faces()
    # extrados = MeshDos.from_vertices_and_faces(vertices, faces)

    i = 0
    for key in intrados.vertices():
        normal_intra = intrados.vertex_attribute(key, 'n')
        normal_extra = extrados.vertex_attribute(key, 'n')
        x, y, z_lb = intrados.vertex_coordinates(key)
        x, y, z_ub = extrados.vertex_coordinates(key)
        if intrados.vertex_attribute(key, '_is_outside'):
            deviation_intra = t
        else:
            deviation_intra = 1/math.sqrt(1/(1 + (normal_intra[0]**2 + normal_intra[1]**2)/normal_intra[2]**2))  # 1/cos(a)
        deviation_extra = 1/math.sqrt(1/(1 + (normal_extra[0]**2 + normal_extra[1]**2)/normal_extra[2]**2))  # 1/cos(a)
        lb_update[i] = z_lb + n * deviation_intra
        ub_update[i] = z_ub - n * deviation_extra
        i += 1

    return ub_update, lb_update


def general_dub_dlb_with_n(ub, lb, n, intrados, extrados, t):

    # ub_update = zeros((len(ub), 1))
    # lb_update = zeros((len(ub), 1))
    dub = zeros((len(ub), 1))
    dlb = zeros((len(ub), 1))

    # vertices, faces = intrados.to_vertices_and_faces()
    # intrados = MeshDos.from_vertices_and_faces(vertices, faces)
    # vertices, faces = extrados.to_vertices_and_faces()
    # extrados = MeshDos.from_vertices_and_faces(vertices, faces)

    # intrados = intrados.copy()
    # extrados = extrados.copy()

    i = 0
    for key in intrados.vertices():
        normal_intra = intrados.vertex_attribute(key, 'n')
        normal_extra = extrados.vertex_attribute(key, 'n')
        x, y, z_lb = intrados.vertex_coordinates(key)
        x, y, z_ub = extrados.vertex_coordinates(key)
        if intrados.vertex_attribute(key, '_is_outside'):
            deviation_intra = t
        else:
            deviation_intra = 1/math.sqrt(1/(1 + (normal_intra[0]**2 + normal_intra[1]**2)/normal_intra[2]**2))  # 1/cos(a)
        deviation_extra = 1/math.sqrt(1/(1 + (normal_extra[0]**2 + normal_extra[1]**2)/normal_extra[2]**2))  # 1/cos(a)
        # lb_update[i] = z_lb + n * deviation_intra
        # ub_update[i] = z_ub - n * deviation_extra
        dlb[i] = + deviation_intra
        dub[i] = - deviation_extra
        i += 1

    # If we must update it before calculating dlb and dub... This is more non linear... Don't actually make sense for the simple thing we are doing....

    # i = 0
    # for key in intrados.vertices():
    #     intrados.vertex_attribute(key, 'z', lb_update[i])
    #     extrados.vertex_attribute(key, 'z', ub_update[i])
    #     i += 1

    # i = 0
    # for key in intrados.vertices():
    #     normal_intra = intrados.vertex_normal(key)
    #     normal_extra = extrados.vertex_normal(key)
    #     x, y, z_lb = intrados.vertex_coordinates(key)
    #     x, y, z_ub = extrados.vertex_coordinates(key)
    #     deviation_intra = 1/math.sqrt(1/(1 + (normal_intra[0]**2 + normal_intra[1]**2)/normal_intra[2]**2))  # 1/cos(a)
    #     deviation_extra = 1/math.sqrt(1/(1 + (normal_extra[0]**2 + normal_extra[1]**2)/normal_extra[2]**2))  # 1/cos(a)
    #     dlb[i] = + deviation_intra
    #     dub[i] = - deviation_extra
    #     i += 1

    return dub, dlb


def _invert_vector(vector):
    return [-1*vector[0], -1*vector[1], -1*vector[2]]


def general_b_update_with_n(b, n, fixed):

    return b_new

def general_db_with_n(b, n, fixed):

    return
