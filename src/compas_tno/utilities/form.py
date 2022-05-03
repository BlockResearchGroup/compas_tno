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
