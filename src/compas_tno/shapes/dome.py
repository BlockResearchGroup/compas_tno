

from scipy.interpolate import interp2d
from numpy import arange
from numpy import zeros
from numpy import array
import math
from math import pi
from math import sin, cos
from compas.datastructures import Mesh
from scipy import interpolate
from scipy import hstack


def set_dome_heighfield(center = [5.0,5.0], radius = 5.0, thk = 0.30, t = 5.0, discretisation=[8, 20]):

    tol = 10e-3
    xc = center[0]
    yc = center[1]
    ri = radius - thk/2
    re = radius + thk/2
    n_radial = discretisation[0]
    n_spikes = discretisation[1]
    r_div = radius/n_radial
    theta = 2*math.pi/n_spikes

    x1d = []
    y1d = []
    zt1d = []
    zi1d = []
    ze1d = []
    faces = []
    faces_i = []
    i = 0
    uv_i = {}

    for nr in range(n_radial+1):
        for nc in range(n_spikes):
            uv_i[(nr, nc)] = i
            xi = xc + nr * r_div * math.cos(theta * nc)
            yi = yc + nr * r_div * math.sin(theta * nc)
            zt2 = radius**2 - (xi - xc)**2 - (yi - yc)**2
            zi2 = ri**2 - (xi - xc)**2 - (yi - yc)**2
            ze2 = re**2 - (xi - xc)**2 - (yi - yc)**2
            ze = math.sqrt(ze2)
            if zt2 < 0.0:
                zt = 0.0
            else:
                zt = math.sqrt(zt2)
            if zi2 < 0.0:
                if zi2 < -tol:
                    zi = -t
                else:
                    zi = 0.0
            else:
                zi = math.sqrt(zi2)

            x1d.append(xi)
            y1d.append(yi)
            zt1d.append(zt)
            zi1d.append(zi)
            ze1d.append(ze)

            if nr < n_radial:
                p1 = (nr, nc)
                p2 = (nr, nc+1)
                p3 = (nr+1, nc)
                p4 = (nr+1, nc+1)
                if nc < n_spikes - 1: # General Case
                    # face = [p1, p2, p4, p3]
                    face = [p3, p4, p2, p1]
                else:
                    p2 = (nr, 0)
                    p4 = (nr+1, 0)
                    # face = [p1, p2, p4, p3]
                    face = [p3, p4, p2, p1]
                faces.append(face)

            i = i+1

    for face in faces:
        face_i = []
        for uv in face:
            u,v = uv
            i = uv_i[(u,v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array([x1d, y1d, zt1d]).transpose()
    xyz_intrados = array([x1d, y1d, zi1d]).transpose()
    xyz_extrados = array([x1d, y1d, ze1d]).transpose()

    intrados = Mesh.from_vertices_and_faces(xyz_intrados, faces_i)
    extrados = Mesh.from_vertices_and_faces(xyz_extrados, faces_i)
    middle = Mesh.from_vertices_and_faces(xyz_middle, faces_i)

    return intrados, extrados, middle


def set_dome_with_spr(center = [5.0,5.0], radius = 5.0, thk = 0.30, theta=[0, math.pi/2], t = 5.0, discretisation=[8, 20]):

    center.append(0.0)
    ri = radius - thk/2
    re = radius + thk/2
    phi_lower = 0
    phi_upper = 2 * math.pi
    phi_length = discretisation[1]
    step_phi = (phi_upper - phi_lower)
    phi_range = [phi_lower + x * (phi_upper - phi_lower) / phi_length for x in range(phi_length + 1)]

    theta_length = discretisation[0]
    theta_lower = theta[0]
    theta_upper = theta[1]
    theta_range = [theta_lower + x * (theta_upper - theta_lower) / theta_length for x in range(theta_length + 1)]

    [xmin, _, _] = geom_dome(center, radius, theta_range[0], phi_range[0])
    [xmax, _, _] = geom_dome(center, radius, theta_range[theta_length], phi_range[phi_length])
    theta_upper_i = (math.pi/2 - math.acos(min(1,(xmax-center[0])/ri)))
    theta_lower_e = (math.pi/2 - math.acos((xmin-center[0])/re))
    print('min/max', xmin, xmax)
    print('lower/lower_i/upper/upper_e', theta_lower/(math.pi/2), theta_lower_e/(math.pi/2), theta_upper/(math.pi/2), theta_upper_i/(math.pi/2))
    theta_range_i = [theta_lower + x * (theta_upper_i - theta_lower) / theta_length for x in range(theta_length + 1)]
    theta_range_e = [theta_lower_e + x * (theta_upper - theta_lower_e) / theta_length for x in range(theta_length + 1)]

    faces = []
    faces_i = []
    count = 0
    uv_i = {}
    xyz_middle = []
    xyz_intrados = []
    xyz_extrados = []

    for i in range(phi_length + 1):
        for j in range(theta_length + 1):
            uv_i[(i, j)] = count
            p1 = (i, j)
            p2 = (i, j+1)
            p3 = (i+1, j)
            p4 = (i+1, j+1)
            face = []

            ptm = geom_dome(center, radius, theta_range[j], phi_range[i])
            pti = geom_dome(center, ri, theta_range_i[j], phi_range[i])
            pte = geom_dome(center, re, theta_range_e[j], phi_range[i])

            xyz_middle.append(ptm)
            xyz_intrados.append(pti)
            xyz_extrados.append(pte)

            if i < phi_length and j < theta_length:
                face = [p3, p4, p2, p1]
                faces.append(face)

            count += 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array(xyz_middle)
    xyz_extrados = array(xyz_extrados)
    xyz_intrados = array(xyz_intrados)

    middle = Mesh.from_vertices_and_faces(xyz_middle, faces_i)
    extrados = Mesh.from_vertices_and_faces(xyz_extrados, faces_i)
    intrados = Mesh.from_vertices_and_faces(xyz_intrados, faces_i)

    return intrados, extrados, middle


def set_dome_polar_coord(center = [5.0,5.0], radius = 5.0, thk = 0.30, theta=[0, math.pi/2], t = 0.0, discretisation=[8, 20]):

    center.append(0.0)
    ri = radius - thk/2
    re = radius + thk/2
    phi_lower = 0
    phi_upper = 2 * math.pi
    phi_length = discretisation[1]
    step_phi = (phi_upper - phi_lower)
    phi_range = [phi_lower + x * (phi_upper - phi_lower) / phi_length for x in range(phi_length + 1)]

    theta_length = discretisation[0]
    theta_lower = theta[0]
    theta_upper = theta[1]
    theta_range = [theta_lower + x * (theta_upper - theta_lower) / theta_length for x in range(theta_length + 1)]

    [xmin, _, _] = geom_dome(center, radius, theta_range[0], phi_range[0])
    [xmax, _, _] = geom_dome(center, radius, theta_range[theta_length], phi_range[phi_length])
    theta_upper_i = (math.pi/2 - math.acos(min(1,(xmax-center[0])/ri)))
    theta_lower_e = (math.pi/2 - math.acos((xmin-center[0])/re))
    print('min/max', xmin, xmax)
    print('lower/lower_i/upper/upper_e', theta_lower/(math.pi/2), theta_lower_e/(math.pi/2), theta_upper/(math.pi/2), theta_upper_i/(math.pi/2))
    theta_range_i = [theta_lower + x * (theta_upper_i - theta_lower) / theta_length for x in range(theta_length + 1)]
    theta_range_e = [theta_lower_e + x * (theta_upper - theta_lower_e) / theta_length for x in range(theta_length + 1)]

    faces = []
    faces_i = []
    count = 0
    uv_i = {}
    xyz_middle = []
    xyz_intrados = []
    xyz_extrados = []

    for i in range(phi_length + 1):
        for j in range(theta_length + 1):
            uv_i[(i, j)] = count
            p1 = (i, j)
            p2 = (i, j+1)
            p3 = (i+1, j)
            p4 = (i+1, j+1)
            face = []

            ptm = geom_dome(center, radius, theta_range[j], phi_range[i])
            pti = geom_dome(center, ri, theta_range_i[j], phi_range[i])
            pte = geom_dome(center, re, theta_range_e[j], phi_range[i])

            xyz_middle.append(ptm)
            xyz_intrados.append(pti)
            xyz_extrados.append(pte)

            if i < phi_length and j < theta_length:
                face = [p3, p4, p2, p1]
                faces.append(face)

            count += 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array(xyz_middle)
    xyz_extrados = array(xyz_extrados)
    xyz_intrados = array(xyz_intrados)

    middle = Mesh.from_vertices_and_faces(xyz_middle, faces_i)
    extrados = Mesh.from_vertices_and_faces(xyz_extrados, faces_i)
    intrados = Mesh.from_vertices_and_faces(xyz_intrados, faces_i)

    return intrados, extrados, middle

def geom_dome(p0, ro, theta, phi):
    x = ro * sin(theta) * cos(phi)
    y = ro * sin(theta) * sin(phi)
    z = ro * cos(theta)
    point = [p0[0] + x, p0[1] +y, p0[2] + z]
    return point
