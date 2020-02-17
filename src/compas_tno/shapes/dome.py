

from scipy.interpolate import interp2d
from numpy import arange
from numpy import zeros
from numpy import array
import math
from compas.datastructures import Mesh
from scipy import interpolate
from scipy import hstack


def set_dome_heighfield(center = [5.0,5.0], radius = 5.0, thk = 0.30, t = 5.0, density=[8, 20]):
    tol = 10e-3
    xc = center[0]
    yc = center[1]
    ri = radius - thk/2
    re = radius + thk/2
    n_radial = density[0]
    n_spikes = density[1]
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
            uv_i[(nr,nc)] = i
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
