from compas_tno.diagrams import FormDiagram
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import distance_point_point_xy


def form_add_lines_support(form, loaded_node, supports):
    """Add direct load paths from a chosen node to the chosen supports

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram to modify
    loaded_node : int
        Key of the node which will be the starting point for the additional lines
    supports : [int]
        List with the keys of the supports that should be linked to the loaded node

    Returns
    -------
    form : FormDiagram
        The form diagram with the modifications

    new_loaded_node : int
        Key for the loaded node at the same position
    """

    text = {}
    new_lines = []
    text[loaded_node] = loaded_node

    xp, yp, _ = form.vertex_coordinates(loaded_node)
    fixed_coords = [form.vertex_coordinates(vertex) for vertex in form.vertices_where({'is_fixed': True})]
    parameters = form.parameters

    lines = form.to_lines()
    support_lines = []
    points = []

    for support in supports:
        xs, ys, _ = form.vertex_coordinates(support)
        support_line = [[xs, ys, 0.0], [xp, yp, 0.0]]
        support_lines.append(support_line)

    for line in lines:
        int_pt1 = intersection_segment_segment_xy(support_lines[0], line)
        int_pt2 = intersection_segment_segment_xy(support_lines[1], line)
        if int_pt1 or int_pt2:
            int_pt = int_pt1 or int_pt2
            points.append(int_pt)
            if int_pt == line[0] or int_pt == line[1]:
                new_lines.append(line)  # necessary?
                pass
            else:
                cutline1 = [line[0], int_pt]
                cutline2 = [int_pt, line[1]]
                new_lines.append(cutline1)
                new_lines.append(cutline2)
        else:
            new_lines.append(line)

    for i in range(len(points) - 1):
        j = i + 1
        if distance_point_point_xy(points[i], points[j]) > 1e-3:
            new_lines.append([points[i], points[j]])

    form = FormDiagram.from_lines(new_lines)

    # find new vertices
    for vertex in form.vertices():
        coord = form.vertex_coordinates(vertex)
        dist = distance_point_point_xy(coord, [xp, yp])
        if dist < 1e-3:
            new_loaded_node = vertex
            continue
        for coord_fix in fixed_coords:
            dist = distance_point_point_xy(coord, coord_fix)
            if dist < 1e-3:
                form.vertex_attribute(vertex, 'is_fixed', True)

    form.parameterss = parameters

    return form, new_loaded_node


def store_inds(form, ind_edges=[]):
    """Store independent edges as an attribute.
    Note: The midpoint of the independent edges is stored.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to store
    ind_edges : list, optional
        List with the independent edges, by default [] in which edges

    Returns
    -------
    None
        Attribute included to the form diagram
    """

    points = []
    if len(ind_edges) == 0:
        for u, v in form.edges_where({'is_ind': True}):
            points.append(Point(*(form.edge_midpoint(u, v)[:2] + [0])))
    else:
        for u, v in ind_edges:
            points.append(Point(*(form.edge_midpoint(u, v)[:2] + [0])))
    form.attributes['indset'] = points

    return


def retrieve_inds(form, indset=[]):
    """Store independent edges as an attribute.
    Note: The midpoint of the independent edges is stored.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to store
    indset : list
        List with the independent edges

    Returns
    -------
    None
        Attribute included to the form diagram
    """

    # WIP
    indset = [a if not isinstance(a, str) else reverse_geometric_key(a) for a in form.attributes['indset']]
    for pt in indset:
        points.append(Point(*(form.edge_midpoint(u, v)[:2] + [0])))
        for u, v in form.edges_where({'_is_edge': True}):
            index = problem.uv_i[(u, v)]
            edgemid = Point(*(form.edge_midpoint(u, v)[:2] + [0]))
            for pt in indset:
                if distance_point_point_xy(edgemid, pt) < tol_old_ind:
                    ind.append(index)
                    break
            if index in ind:
                indset.remove(pt)
    form.attributes['indset'] = points

    return


def retrieve_form_polylines(form, singularities=[]):
    """Store independent edges as an attribute.
    Note: The midpoint of the independent edges is stored.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to store
    singularities : list
        List with index of singularities

    Returns
    -------
    polylines : list[list]
        List with the edges forming polylines
    """

    polylines = []
    edges_visited = []
    edges = list(form.edges_where({'_is_edge': True}))

    for edge in edges:
        if edge in edges_visited:
            continue
        edge

    return


def mesh_remove_two_valent_nodes(mesh, cls=Mesh, delete_boundary_face=True):
    """Remove 2 valent nodes from a mesh or form diagram

    Parameters
    ----------
    mesh : Mesh
        Mesh to perform operations
    cls : class, optional
        Class to return the mesh, by default Mesh

    Returns
    -------
    mesh
        Mesh without nodes
    """

    uv_i = {(u, v): index for index, (u, v) in enumerate(mesh.edges())}
    delete_indices = []
    newedges = []
    for key in mesh.vertices():
        deg = mesh.vertex_degree(key)
        if deg == 2:
            pts = mesh.vertex_neighbors(key)
            ed = mesh.vertex_edges(key)
            delete_indices.append(uv_i[ed[0]])
            delete_indices.append(uv_i[ed[1]])
            newedges.append([mesh.vertex_coordinates(a) for a in pts])

    lines = mesh.to_lines()
    new_lines = []

    for i, line in enumerate(lines):
        if i in delete_indices:
            continue
        new_lines.append(line)

    for edge in newedges:
        new_lines.append(edge)

    mesh = cls.from_lines(new_lines, delete_boundary_face=delete_boundary_face)

    return mesh


def form_parabolic_slide(delta, y0, y1, form):
    """Modify the form diagram applying a parabolic displacement profile to the nodes

    Parameters
    ----------
    delta : float
        Maximum distance applied to the nodes
    y0 : float
        Start ordinate to apply the parabolic pattern
    y1 : float
        End ordinate to apply the parabolic pattern
    form : FormDiagram
        The Form Diagram

    Returns
    -------
    form : FormDiagram
        The form diagram with the modifications

    """

    yc = ((y1 - y0)/2)
    for vertex in form.vertices():
        x, y, _ = form.vertex_coordinates(vertex)
        dy = min(y - y0, y1 - y)
        if abs(dy) > 1e-3:
            dx = delta * (1 - ((dy - yc)/yc)**2)
            form.vertex_attribute(vertex, 'x', x + dx)

    return form


def move_pattern_to_origin(mesh, corners=[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]], fix_corners=True):
    """Modify the form diagram to coincide the lower corner of its bounding box with the corners provided. Note: translation and scaling are involved.

    Parameters
    ----------
    mesh : Mesh
        The form diagram adopted
    corners : [pts]
        A list with four points for the corner

    Returns
    -------
    trans : list
        List of the translation factors
    factors : list
        List of the scaling factors
    """

    xmax_, ymax_ = max([a[0] for a in corners]), max([a[1] for a in corners])
    tol = 10e-3

    trans = []
    factors = []

    bbox = mesh.bounding_box_xy()
    xmin, ymin = min([m[0] for m in bbox]), min([m[1] for m in bbox])
    # print('bbox beginning', bbox)

    if abs(xmin) > tol or abs(ymin) > tol:
        dx, dy = -xmin, -ymin
        trans = [dx, dy, 0.0]
        translation = Translation.from_vector(trans)
        mesh.transform(translation)

    bbox = mesh.bounding_box_xy()
    xmax, ymax = max([m[0] for m in bbox]), max([m[1] for m in bbox])
    # print('bbox 2', bbox)
    # print(xmax, ymax)

    if abs(ymax_ - ymax) > tol or abs(xmax_ - xmax) > tol:
        factors = [xmax_/xmax, ymax_/ymax, 0.0]
        scale = Scale.from_factors(factors)
        mesh.transform(scale)

    if fix_corners:
        for key in mesh.vertices():
            pt = mesh.vertex_coordinates(key)
            for corner in corners:
                dist = distance_point_point_xy(pt, corner)
                if dist < 1e-3:
                    mesh.vertex_attribute(key, 'is_fixed', True)
                    break

    return trans, factors


def slide_diagram(form, delta=0.5, y0=0.0, y1=10.0):
    """Apply a parabolic sliding to the nodes towards +x direction. Sliding profile is defined upon pattern's height (y-coordinate)

    Parameters
    ----------
    form : The form diagram of the problem
        _description_
    delta : float, optional
        Delta distance to apply to the nodes, by default 0.5
    y0 : float, optional
        Start of parabolic profile, by default 0.0
    y1 : float, optional
        End of parabolic profile, by default 10.0
    """

    yc = ((y1 - y0)/2)
    for vertex in form.vertices_where({'is_fixed': False}):
        x, y, _ = form.vertex_coordinates(vertex)
        dy = min(y - y0, y1 - y)
        if abs(dy) > 1e-3:
            dx = delta * (1 - ((dy - yc)/yc)**2)
            form.vertex_attribute(vertex, 'x', x + dx)


def displacement_map_parabola(form, y0=0.0, y1=10.0):
    """Create the displacement map for a 1D parabola sliding of the structural pattern.

    Parameters
    ----------
    form : The form diagram of the problem
        _description_
    delta : float, optional
        Delta distance to apply to the nodes, by default 0.5
    y0 : float, optional
        Start of parabolic profile, by default 0.0
    y1 : float, optional
        End of parabolic profile, by default 10.0

    Returns
    -------
    dX: array [n, 2]
        Displacement field in plan.
    """

    n = form.number_of_vertices()
    k = - 4/(y1 - y0)**2
    dX = zeros((n, 2))

    for i, vertex in enumerate(form.vertices()):
        if form.vertex_attribute(vertex, 'is_fixed'):
            continue
        x, y, _ = form.vertex_coordinates(vertex)
        dyi = 0.0
        dxi = k * (y - y0) * (y - y1)
        dX[i] = [dxi, dyi]

    return dX


def displacement_map_4parabolas(form, y0=0.0, y1=10.0, x0=0.0, x1=10.0, tol=0.1):
    """Create the displacement map based on a paraboloid sliding the structural pattern.

    Parameters
    ----------
    form : The form diagram of the problem
        _description_
    delta : float, optional
        Delta distance to apply to the nodes, by default 0.5
    y0 : float, optional
        Start of parabolic profile, by default 0.0
    y1 : float, optional
        End of parabolic profile, by default 10.0

    Returns
    -------
    dX: array [n, 2]
        Displacement field in plan.
    """

    n = form.number_of_vertices()
    kx = - 4/(y1 - y0) ** 2
    ky = - 4/(x1 - x0) ** 2
    dX = zeros((n, 2))

    for i, vertex in enumerate(form.vertices()):
        if form.vertex_attribute(vertex, 'is_fixed'):
            continue
        x, y, _ = form.vertex_coordinates(vertex)
        dxi = 0.0
        dyi = 0.0

        if y <= y0 + (y1 - y0)/(x1 - x0) * (x - x0) - tol and y >= y1 - (y1 - y0)/(x1 - x0) * (x - x0) + tol:  # Q1
            dxi = - kx * (y - y0) * (y - y1)
        elif y >= y0 + (y1 - y0)/(x1 - x0) * (x - x0) + tol and y >= y1 - (y1 - y0)/(x1 - x0) * (x - x0) + tol:  # Q3
            dyi = - ky * (x - x0) * (x - x1)
        elif y >= y0 + (y1 - y0)/(x1 - x0) * (x - x0) + tol and y <= y1 - (y1 - y0)/(x1 - x0) * (x - x0) - tol:  # Q2
            dxi = + kx * (y - y0) * (y - y1)
        elif y <= y0 + (y1 - y0)/(x1 - x0) * (x - x0) - tol and y <= y1 - (y1 - y0)/(x1 - x0) * (x - x0) - tol:  # Q4
            dyi = ky * (x - x0) * (x - x1)
        else:  # diagonals
            pass

        dX[i] = [dxi, dyi]

    return dX


def displacement_map_paraboloid(form, xc=5.0, yc=5.0, radius=5.0):
    """Create the displacement map based on a paraboloid sliding the structural pattern.

    Parameters
    ----------
    form : The form diagram of the problem
        _description_
    delta : float, optional
        Delta distance to apply to the nodes, by default 0.5
    y0 : float, optional
        Start of parabolic profile, by default 0.0
    y1 : float, optional
        End of parabolic profile, by default 10.0

    Returns
    -------
    dX: array [n, 2]
        Displacement field in plan.
    """

    n = form.number_of_vertices()
    k = 1/radius**2
    dX = zeros((n, 2))

    for i, vertex in enumerate(form.vertices()):
        if form.vertex_attribute(vertex, 'is_fixed'):
            continue
        x, y, _ = form.vertex_coordinates(vertex)
        dxi = k * (x - xc) ** 2 * 2 * (xc - x)/ radius
        dyi = k * (y - yc) ** 2 * 2 * (yc - y)/ radius

        dX[i] = [dxi, dyi]

    return dX
