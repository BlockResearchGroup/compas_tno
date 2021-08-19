import math
from compas.geometry import subtract_vectors
from compas.geometry import normalize_vector
from compas.geometry import angle_vectors
from compas.geometry import scale_vector
from compas.geometry import Vector
from compas.geometry import rotate_points
from compas_assembly.datastructures import Block
from compas_rhino.artists import MeshArtist
from compas_viewers.multimeshviewer import MultiMeshViewer


import os
from math import pi
from compas.geometry import Rotation
from compas_assembly.datastructures import Assembly
from compas_assembly.datastructures import assembly_interfaces_numpy
from compas_assembly.plotter import AssemblyPlotter


__all__ = [
    'semicirc_arch_height_span',
]


def semicirc_arch_height_span(h, s, d, t, n_blocks):
    """Create a semicircular arch from the dimension of its height and span.

        Parameters
        ----------
        h : float
            Rise (Height) of the arch in meter measured from the impost to the intrados
            of the keystone.
        s : float
            Dimension of the span in meter measured at the impost level.
        d : float
            Depth in meter of the arch perpendicular to the front view plane.
        t : float
            Thickness in meter of the arch.
        n_blocks: int
            number of voussoirs


        Returns
        -------
        list of float
            The XYZ coordinates of the closest point.
"""

    if h > s / 2:
        print("Not a semicircular arch")
        exit()

    radius = h / 2 + (s**2 / (8 * h))
    print('radius =', radius)

    pt0 = (0.0, 0.0, 0.0)
    pt1 = (0.0, 0.0, h)
    pt2 = (pt0[0] - s / 2, 0.0, 0.0)

    radius_pt = (pt0[0], pt0[1], pt0[2] - (radius - pt1[2]))
    vect = normalize_vector(subtract_vectors(pt2, radius_pt))
    angle_vect = angle_vectors(vect, scale_vector(Vector.Xaxis(), -1.0))
    springing_angle = math.radians(90) - angle_vect
    print('springing_angle =', math.degrees(springing_angle))
    tot_angle = math.radians(180) - 2 * angle_vect

    a = pt1
    b = (pt1[0], pt1[1] + d, pt1[2])
    c = (pt1[0], pt1[1] + d, pt1[2] + t)
    d = (pt1[0], pt1[1], pt1[2] + t)

    an = tot_angle / n_blocks

    voussoirs = []
    pts = rotate_points([a, b, c, d], springing_angle, scale_vector(Vector.Yaxis(), -1.0), radius_pt)
    points = []
    faces = [[0, 1, 2, 3], [0, 3, 7, 4], [0, 4, 5, 1], [5, 6, 2, 1], [6, 7, 3, 2], [5, 4, 7, 6]]

    for n in range(n_blocks + 1):
        pts1 = rotate_points([pts[0], pts[1], pts[2], pts[3]], an * n, Vector.Yaxis(), radius_pt)
        points.append(pts1)
    for n in range(n_blocks):
        m = Block.from_vertices_and_faces((points[n] + points[n + 1]), faces)
        voussoirs.append(m)

    dist = 0.25
    for k in [0, -1]:
        arch_face = points[k]
        block_base = [point[:] for point in points[k]]
        for pt in block_base:
            pt[2] = - dist
        m = Block.from_vertices_and_faces((arch_face + block_base), faces)
        voussoirs.append(m)

    return voussoirs

if __name__ == "__main__":

    h = 2
    s = 10
    d = 0.7
    t = 1
    n_blocks = 10

    voussoirs = semicirc_arch_height_span(h, s, d, t, n_blocks)

    viewer = MultiMeshViewer()
    viewer.meshes = voussoirs
    viewer.show()

    # R = Rotation.from_axis_and_angle([1.0, 0, 0], -pi / 2)
    # assembly.transform(R)

    # plotter = AssemblyPlotter(assembly, figsize=(16, 10), tight=True)

    # plotter.draw_nodes(radius=0.02, facecolor={key: "#ff0000" for key in assembly.nodes_where({'is_support': True})})
    # plotter.draw_edges()
    # plotter.draw_blocks(edgecolor='#444444', edgewidth=0.5)
    # plotter.show()

