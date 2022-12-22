from compas.geometry import distance_point_point_xy
from compas.geometry import distance_point_line_xy
from compas.geometry import closest_point_on_line_xy

from numpy import zeros

import matplotlib.pyplot as plt


__all__ = [
    'apply_radial_symmetry',
    'apply_symmetry_from_axis',
    'find_sym_axis_in_rect_patterns',
    'build_symmetry_matrix',
    'build_symmetry_transformation',
    'build_vertex_symmetry_transformation',
    'build_symmetry_matrix_supports'
]


def apply_radial_symmetry(form, center=[5.0, 5.0, 0.0], correct_loads=True):
    """Apply a radial symmetry based on a center points. Applicable for dome circulat patterns.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    center : [float, float, float], optional
        The coordinates of the center of the pattern, by default [5.0, 5.0, 0.0]
    correct_loads : bool, optional
        Whether or not the loads should be corrected in the nodes for perfect symmetry, by default True

    Returns
    -------
    None
        The FormDiagram is modified in place.
    """

    form.edges_attribute('sym_key', None)
    form.vertices_attribute('sym_key', None)

    dist_checked = []
    dist_dict = {}

    # Symmetry on the independent edges

    for u, v in form.edges_where({'_is_edge': True}):
        midpoint = form.edge_midpoint(u, v)
        dist = round(distance_point_point_xy(center, midpoint), 10)
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

    for key in form.vertices():
        point = form.vertex_coordinates(key)
        dist = round(distance_point_point_xy(center, point), 10)
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
        i += 1

    return


def apply_symmetry_from_axis(form, list_axis_symmetry=[], correct_loads=True, tol=0.01):
    """ Apply a symmetry based on a series of axis of symmetry. Applicable for rectangular patterns.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    list_axis_symmetry : list, optional
        The list of the axis of symmetry, by default []
    correct_loads : bool, optional
        Whether or not the loads should be corrected in the nodes for perfect symmetry, by default True
    tol : float, optional
        Tollerance to assume symmetry, by default 0.001

    Returns
    -------
    None
        The FormDiagram is modified in place.

    """

    form.edges_attribute('sym_dict', None)
    form.vertices_attribute('sym_dict', None)
    form.edges_attribute('sym_key', None)
    form.vertices_attribute('sym_key', None)

    i = 0
    for axis_symmetry in list_axis_symmetry:
        axis_str = str(axis_symmetry)
        dist_checked = []
        dist_dict = {}
        for u, v in form.edges_where({'_is_edge': True}):
            midpoint = form.edge_midpoint(u, v)
            dist_line = round(distance_point_line_xy(midpoint, axis_symmetry), 10)  # change this logic to something that takes into account tol
            closest_pt = closest_point_on_line_xy(midpoint, axis_symmetry)  # geometric_key(..) used before
            # dist = [closest_pt, dist_line]
            dist = [dist_line, closest_pt]
            dist_dict[(u, v)] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
        for dist in dist_checked:
            for u, v in dist_dict:
                # if dist_dict[(u, v)] == dist:
                if abs(dist[0] - dist_dict[(u, v)][0]) < tol and distance_point_point_xy(dist[1], dist_dict[(u, v)][1]) < tol:
                    form.edge_attribute((u, v), 'sym_key', i)
                    if not form.edge_attribute((u, v), 'sym_dict'):
                        dic = {}
                    else:
                        dic = form.edge_attribute((u, v), 'sym_dict')
                    dic[axis_str] = i
                    form.edge_attribute((u, v), 'sym_dict', dic)
            i += 1

    if len(list_axis_symmetry) > 1:
        checked_edge = []
        groups_of_edges = []
        for edge in form.edges_where({'_is_edge': True}):
            if edge in checked_edge:
                continue
            group = [edge]
            checked_edge.append(edge)
            values = list(form.edge_attribute(edge, 'sym_dict').values())
            for edge2 in form.edges_where({'_is_edge': True}):
                if edge2 in checked_edge:
                    continue
                values2 = list(form.edge_attribute(edge2, 'sym_dict').values())
                if len(set(values2) & set(values)):
                    group.append(edge2)
                    checked_edge.append(edge2)
                    values = list(set(values + values2))
            groups_of_edges.append(group)

        i = 0
        for group in groups_of_edges:
            form.edges_attribute('sym_key', value=i, keys=group)
            i += 1

    i = 0
    for axis_symmetry in list_axis_symmetry:
        axis_str = str(axis_symmetry)
        dist_checked = []
        dist_dict = {}
        for key in form.vertices():
            point = form.vertex_coordinates(key)
            dist_line = round(distance_point_line_xy(point, axis_symmetry), 10)
            closest_pt = closest_point_on_line_xy(point, axis_symmetry)  # geometric_key( .. )
            # dist = [closest_pt, dist_line]
            dist = [dist_line, closest_pt]
            dist_dict[key] = dist
            if dist not in dist_checked:
                dist_checked.append(dist)
        for dist in dist_checked:
            for key in dist_dict:
                # if dist_dict[key] == dist:
                if abs(dist[0] - dist_dict[key][0]) < tol and distance_point_point_xy(dist[1], dist_dict[key][1]) < tol:
                    form.vertex_attribute(key, 'sym_key', i)
                    if not form.vertex_attribute(key, 'sym_dict'):
                        dic = {}
                    else:
                        dic = form.vertex_attribute(key, 'sym_dict')
                    dic[axis_str] = i
                    form.vertex_attribute(key, 'sym_dict', dic)
            i += 1

    if len(list_axis_symmetry) > 1:
        checked_vertex = []
        groups_of_vertices = []
        for vertex in form.vertices():
            if vertex in checked_vertex:
                continue
            group = [vertex]
            checked_vertex.append(vertex)
            values = list(form.vertex_attribute(vertex, 'sym_dict').values())
            for vertex2 in form.vertices():
                if vertex2 in checked_vertex:
                    continue
                values2 = list(form.vertex_attribute(vertex2, 'sym_dict').values())
                if len(set(values2) & set(values)):
                    group.append(vertex2)
                    checked_vertex.append(vertex2)
                    values = list(set(values + values2))
            groups_of_vertices.append(group)

        i = 0
        for group in groups_of_vertices:
            form.vertices_attribute('sym_key', value=i, keys=group)
            i += 1
            if correct_loads:
                pz = form.vertex_attribute(group[0], 'pz')
                ub = form.vertex_attribute(group[0], 'ub')
                lb = form.vertex_attribute(group[0], 'lb')
                s = form.vertex_attribute(group[0], 'target')
                form.vertices_attribute('pz', value=pz, keys=group)
                form.vertices_attribute('ub', value=ub, keys=group)
                form.vertices_attribute('lb', value=lb, keys=group)
                form.vertices_attribute('target', value=s, keys=group)

    return


def find_sym_axis_in_rect_patterns(data_form):
    """ Find the axis of symmetry in rectangular patterns.

    Parameters
    ----------
    data_form : dict
        A dictionary with the parameters used to create the FormDiagram.

    Returns
    -------
    lines: list
        A list of the symmetry lines.

    """

    lines = []
    xy_span = data_form.get('xy_span', None)
    if xy_span:
        [[x0, x1], [y0, y1]] = xy_span
        xc = (x0 + x1)/2
        yc = (y0 + y1)/2
        hor_line = [[x0, yc], [x1, yc]]
        ver_line = [[xc, y0], [xc, y1]]
        diag_line = [[x0, y0], [x1, y1]]
        lines = [hor_line, ver_line, diag_line]
    else:
        raise NotImplementedError

    return lines


def build_symmetry_matrix(form, printout=False):
    """ Build a symmetry matrix such as Asym * q = 0, with Asym shape (m - k; m)

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    printout : bool, optional
        Whether or not display messages are printed, by default True

    Returns
    -------
    Asym: array (m - k x m)
        The symmetry matrix.
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
    """Build a symmetry matrix Esym (m, k) such as q = Esym * qsym.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    printout : bool, optional
        Whether or not display messages are printed, by default True

    Returns
    -------
    Esym: array (m x k)
        The symmetry matrix.

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
    """Build a symmetry matrix Evsym (n, k) such as z = Evsym * zb.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    printout : bool, optional
        Whether or not display messages are printed, by default True

    Returns
    -------
    Evsym: array (n x k)
        The symmetry matrix.
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
    """Build a symmetry matrix to the supports.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram of the problem
    printout : bool, optional
        Whether or not display messages are printed, by default True

    Returns
    -------
    Evsym: array (n x k)
        The symmetry matrix.
    """

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
