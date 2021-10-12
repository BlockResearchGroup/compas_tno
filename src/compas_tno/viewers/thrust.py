from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import Point

from compas_tno.shapes import Shape

from compas_view2 import app
from compas_view2.shapes import Arrow

from compas_tno.algorithms import reactions


__all__ = [
    'view_thrust',
    'view_solution'
]


def view_thrust(form, settings_form=None, cracks=True, show_reactions=True, thickness=True, outside=True):
    """ Viewer showing the thrust network.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    cracks : bool
        Option to show the cracks as points.
        The default value is ``True`` in which case points touching intrados and exterados are highlighted
    show_reactions : bool
        Whether or not the reaction forces are displayed.
        The default value is ``True``.
    thickness : bool
        Whether or not scale the thickness of the lines in the thrust network to the forces carried.
        The default value is ``True``.
    outside : bool
        Whether or not the nodes outside the envelope are highlighted in black.
        The default value is ``True``.

    Returns
    ----------
    obj
        Viewer object.

    """

    reactions(form)

    viewer = app.App()

    # modify default settings - should be improved
    viewer.view.camera.target = [0, 0, 0]
    viewer.view.camera.distance = 40
    viewer.view.camera.rx = -45
    viewer.view.camera.rz = 45
    viewer.view.camera.fov = 40

    intrad = 1
    extrad = 1
    out = 1

    for key in form.vertices():
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        x, y, z = form.vertex_coordinates(key)
        if cracks:
            if abs(ub - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Extrados (%s)" % extrad, color=(0, 0.5, 0), size=15.0)  # green extrados
                extrad += 1
            elif abs(lb - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Intrados (%s)" % intrad, color=(0, 0, 1), size=15.0)   # blue intrados
                intrad += 1
        if outside:
            if z > ub:
                viewer.add(Point(x, y, z), name="Outside - Intra (%s)" % out, color=(0, 0, 0), size=15.0)   # black outside
                out += 1
            elif z < lb:
                viewer.add(Point(x, y, z), name="Outside - Extra (%s)" % out, color=(0, 0, 0), size=15.0)   # black outside
                out += 1

    forces = [form.edge_attribute((u, v), 'q') * form.edge_length(u, v) for u, v in form.edges_where({'_is_edge': True})]
    fmax = max(abs(max(forces)), abs(min(forces)))
    max_thick = 10.0

    for u, v in form.edges_where({'_is_edge': True}):
        Xu = form.vertex_coordinates(u)
        Xv = form.vertex_coordinates(v)
        line = Line(Xu, Xv)
        if not thickness:
            viewer.add(line, name=str((u, v)), linewidth=max_thick, color=(255, 0, 0))
            continue
        q = form.edge_attribute((u, v), 'q')
        length = form.edge_length(u, v)
        force = abs(q*length)
        thk = force/fmax * max_thick
        if force > 10e-4:
            viewer.add(line, name=str((u, v)), linewidth=thk, color=(255, 0, 0))

    if show_reactions:
        reaction_scale = 0.005

        for key in form.vertices_where({'is_fixed': True}):
            x, y, z = form.vertex_coordinates(key)
            rx = form.vertex_attribute(key, '_rx') * reaction_scale
            ry = form.vertex_attribute(key, '_ry') * reaction_scale
            rz = form.vertex_attribute(key, '_rz') * reaction_scale
            arrow = Arrow([x, y, z], [-rx, -ry, -rz], head_width=0.1, body_width=0.01)
            viewer.add(arrow, color=(0.8, 0.8, 0.8))

    return viewer


def view_solution(form, shape=None, cracks=True, show_reactions=True, thickness=True, outside=True):
    """ Viewer showing the solution: Thrust Network and Masonry Shape.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    shape : Shape
        The Masonry Shape to plot.
        The default value is ``None`` in which case the Shape is derived from the form's attributes
    cracks : bool
        Option to show the cracks as points.
        The default value is ``True`` in which case points touching intrados and exterados are highlighted
    show_reactions : bool
        Whether or not the reaction forces are displayed.
        The default value is ``True``.
    thickness : bool
        Whether or not scale the thickness of the lines in the thrust network to the forces carried.
        The default value is ``True``.
    outside : bool
        Whether or not the nodes outside the envelope are highlighted in black.
        The default value is ``True``.

    Returns
    ----------
    obj
        Viewer object.

    """

    reactions(form)

    if not shape:
        shape = Shape.from_formdiagram_and_attributes(form)

    viewer = app.App()

    # modify default settings - should be improved
    viewer.view.camera.target = [0, 0, 0]
    viewer.view.camera.distance = 40
    viewer.view.camera.rx = -45
    viewer.view.camera.rz = 45
    viewer.view.camera.fov = 40

    intrad = 1
    extrad = 1
    out = 1

    for key in form.vertices():
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        x, y, z = form.vertex_coordinates(key)
        if cracks:
            if abs(ub - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Extrados (%s)" % extrad, color=(0, 0.5, 0), size=15.0)  # green extrados
                extrad += 1
            elif abs(lb - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Intrados (%s)" % intrad, color=(0, 0, 1), size=15.0)   # blue intrados
                intrad += 1
        if outside:
            if z > ub:
                viewer.add(Point(x, y, z), name="Outside - Intra (%s)" % out, color=(0, 0, 0), size=15.0)   # black outside
                out += 1
            elif z < lb:
                viewer.add(Point(x, y, z), name="Outside - Extra (%s)" % out, color=(0, 0, 0), size=15.0)   # black outside
                out += 1

    forces = [form.edge_attribute((u, v), 'q') * form.edge_length(u, v) for u, v in form.edges_where({'_is_edge': True})]
    fmax = max(abs(max(forces)), abs(min(forces)))
    max_thick = 10.0

    for u, v in form.edges_where({'_is_edge': True}):
        Xu = form.vertex_coordinates(u)
        Xv = form.vertex_coordinates(v)
        line = Line(Xu, Xv)
        if not thickness:
            viewer.add(line, name=str((u, v)), linewidth=max_thick, color=(255, 0, 0))
            continue
        q = form.edge_attribute((u, v), 'q')
        length = form.edge_length(u, v)
        force = abs(q*length)
        thk = force/fmax * max_thick
        if force > 10e-4:
            viewer.add(line, name=str((u, v)), linewidth=thk, color=(255, 0, 0))

    if show_reactions:
        reaction_scale = 0.005

        for key in form.vertices_where({'is_fixed': True}):
            x, y, z = form.vertex_coordinates(key)
            rx = form.vertex_attribute(key, '_rx') * reaction_scale
            ry = form.vertex_attribute(key, '_ry') * reaction_scale
            rz = form.vertex_attribute(key, '_rz') * reaction_scale
            arrow = Arrow([x, y, z], [-rx, -ry, -rz], head_width=0.1, body_width=0.01)
            viewer.add(arrow, color=(0.8, 0.8, 0.8))

    vertices_intra, faces_intra = shape.intrados.to_vertices_and_faces()
    mesh_intra = Mesh.from_vertices_and_faces(vertices_intra, faces_intra)

    vertices_extra, faces_extra = shape.extrados.to_vertices_and_faces()
    mesh_extra = Mesh.from_vertices_and_faces(vertices_extra, faces_extra)

    viewer.add(mesh_intra, name="Intrados", show_edges=False, opacity=0.5, color=(0.5, 0.5, 0.5))
    viewer.add(mesh_extra, name="Extrados", show_edges=False, opacity=0.5, color=(0.5, 0.5, 0.5))

    return viewer
