from compas.utilities import geometric_key

from compas.geometry import distance_point_point_xy
from compas.geometry import distance_point_line_xy
from compas.geometry import closest_point_on_line_xy

from numpy import zeros

import matplotlib.pyplot as plt


__all__ = [
    'apply_symmetry',
    'build_symmetry_matrix',
    'build_symmetry_transformation',
    'build_vertex_symmetry_transformation',
    'build_symmetry_matrix_supports'
]


def apply_symmetry(form, center=[5.0, 5.0, 0.0], axis_symmetry=None, correct_loads=True):

    form.edges_attribute('sym_key', None)
    form.vertices_attribute('sym_key', None)
    dist_checked = []
    dist_dict = {}

    # Symmetry on the independent edges

    if not axis_symmetry:
        for u, v in form.edges_where({'_is_edge': True}):
            midpoint = form.edge_midpoint(u, v)
            dist = round(distance_point_point_xy(center, midpoint), 10)
            dist_dict[(u, v)] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
    else:
        for u, v in form.edges_where({'_is_edge': True}):
            midpoint = form.edge_midpoint(u, v)
            dist_line = round(distance_point_line_xy(midpoint, axis_symmetry), 10)
            closest_pt = geometric_key(closest_point_on_line_xy(midpoint, axis_symmetry))
            dist = [closest_pt, dist_line]
            dist_dict[(u, v)] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
    i = 0
    for dist in dist_checked:
        for u, v in dist_dict:
            if dist_dict[(u, v)] == dist:
                form.edge_attribute((u, v), 'sym_key', i)

        i += 1

    # Symmetry on the Support's position

    dist_checked = []
    dist_dict = {}

    if not axis_symmetry:
        for key in form.vertices():
            point = form.vertex_coordinates(key)
            dist = round(distance_point_point_xy(center, point), 10)
            dist_dict[key] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
    else:
        for key in form.vertices():
            point = form.vertex_coordinates(key)
            dist_line = round(distance_point_line_xy(point, axis_symmetry), 10)
            closest_pt = geometric_key(closest_point_on_line_xy(point, axis_symmetry))
            dist = [closest_pt, dist_line]
            dist_dict[key] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
    i = 0
    for dist in dist_checked:
        pass_first = False
        pz = None
        ub = None
        lb = None
        for key in dist_dict:
            if dist_dict[key] == dist:
                form.vertex_attribute(key, 'sym_key', i)
                if correct_loads:
                    if not pass_first:
                        pass_first = True
                        pz = form.vertex_attribute(key, 'pz')
                        ub = form.vertex_attribute(key, 'ub')
                        lb = form.vertex_attribute(key, 'lb')
                        s = form.vertex_attribute(key, 'target')
                    form.vertex_attribute(key, 'pz', pz)
                    form.vertex_attribute(key, 'ub', ub)
                    form.vertex_attribute(key, 'lb', lb)
                    form.vertex_attribute(key, 'target', s)
                    # print(pz, ub, lb, s)
        i += 1

    return


def build_symmetry_matrix(form, printout=False):
    """
    Build a symmetry matrix such as Asym * q = 0, with Asym shape (m - k; m)
    """

    m = len(list(form.edges_where({'_is_edge': True})))
    k_unique = form.number_of_sym_edges(printout=printout)
    Asym = zeros((m - k_unique, m))
    uv_i = form.uv_index()

    line = 0
    for id_sym in range(k_unique):
        i = 0
        for u, v in form.edges():
            if form.edge_attribute((u, v), 'sym_key') == id_sym:
                index = uv_i[(u, v)]
                if i == 0:
                    index0 = index
                else:
                    Asym[line, index0] = 1
                    Asym[line, index] = -1
                    line += 1
                i += 1
    if printout:
        plt.matshow(Asym)
        plt.colorbar()
        plt.show()

    return Asym


def build_symmetry_transformation(form, printout=False):
    r"""
    Build a symmetry matrix Esym (m, k) such as q = Esym * qsym.
    """

    m = len(list(form.edges_where({'_is_edge': True})))
    k_unique = form.number_of_sym_edges(printout=printout)
    Esym = zeros((m, k_unique))
    uv_i = form.uv_index()

    for id_sym in range(k_unique):
        Ei = zeros((m, 1))
        for u, v in form.edges():
            if form.edge_attribute((u, v), 'sym_key') == id_sym:
                index = uv_i[(u, v)]
                Ei[index] = 1.0
        Esym[:, id_sym] = Ei.flatten()

    if printout:
        plt.matshow(Esym)
        plt.colorbar()
        plt.show()

    return Esym


def build_vertex_symmetry_transformation(form, printout=False):
    r"""
    Build a symmetry matrix Evsym (n, k) such as z = Evsym * z_.
    """

    n = form.number_of_edges()
    k_unique = form.number_of_sym_vertices(printout=printout)
    Evsym = zeros((n, k_unique))
    k_i = form.key_index()

    for id_sym in range(k_unique):
        Ei = zeros((n, 1))
        for key in form.vertices():
            if form.vertex_attribute(key, 'sym_key') == id_sym:
                index = k_i[key]
                Ei[index] = 1.0
        Evsym[:, id_sym] = Ei.flatten()

    if printout:
        plt.matshow(Evsym)
        plt.colorbar()
        plt.show()

    return Evsym


def build_symmetry_matrix_supports(form, printout=False):

    n = form.number_of_supports()
    k_unique = form.number_of_sym_supports(printout=printout)
    Asym = zeros((n - k_unique, n))
    key_i_sup = {}

    i = 0
    for key in form.vertices_where({'is_fixed': True}):
        key_i_sup[key] = i
        i += 1

    line = 0
    for id_sym in range(k_unique):
        i = 0
        for key in form.vertices_where({'is_fixed': True}):
            if form.vertex_attribute(key, 'sym_key') == id_sym:
                index = key_i_sup[key]
                if i == 0:
                    index0 = index
                else:
                    Asym[line, index0] = 1
                    Asym[line, index] = -1
                    line += 1
                i += 1
    if printout:
        plt.matshow(Asym)
        plt.colorbar()
        plt.show()

    return Asym
