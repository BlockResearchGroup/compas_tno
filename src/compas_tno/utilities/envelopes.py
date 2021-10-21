from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.dome import dome_zt_update

from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_middle_update

from compas_tno.shapes.pointed_crossvault import pointed_vault_ub_lb_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update

from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.pointed_arch import pointed_arch_ub_lb_update

from compas_tno.utilities.interpolation import get_shape_ub_pattern
from compas_tno.utilities.interpolation import get_shape_lb_pattern
from compas_tno.utilities.interpolation import get_shape_middle_pattern


__all__ = [
    'apply_envelope_from_shape',
    'apply_envelope_on_xy',
    'apply_bounds_on_q',
    'project_mesh_to_middle',
    'apply_envelope_from_shape_proxy',
]


def apply_envelope_from_shape_proxy(formdata, shapedata):
    """ Normalise the color in RGB dividing each component by 255

    Parameters
    ----------
    rgb : tuple
        Tuple with RBG colors.

    Returns
    ----------
    norm_rgb : tuple
        The normalised color.
    """

    from compas_tno.diagrams import FormDiagram
    from compas_tno.shapes import Shape

    form = FormDiagram.from_data(formdata)
    shape = Shape.from_data(shapedata)

    apply_envelope_from_shape(form, shape)

    return form.to_data()


def apply_envelope_from_shape(form, shape):
    """ Apply an envelope (intrados and extrados) to the FormDiagram based on the input shape.

    Parameters
    ----------
    form : ::FormDiagram::
        The input FormDiagram.
    shape : ::Shape::
        The input Shape with intrados and extrardos.

    Returns
    ----------
    None
        The formdiagram is updated in place.
    """

    x = form.vertices_attribute('x')  # check if array is necessary here
    y = form.vertices_attribute('y')

    if shape.datashape['type'] == 'dome':
        zub, zlb = dome_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['center'], shape.datashape['radius'])
    elif shape.datashape['type'] == 'crossvault':
        zub, zlb = crossvault_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        zub, zlb = pointed_vault_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['xy_span'], hc=shape.datashape['hc'],
                                              he=shape.datashape['he'], hm=shape.datashape['hm'])
    elif shape.datashape['type'] == 'arch':
        zub, zlb = arch_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], H=shape.datashape['H'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'pointed_arch':
        zub, zlb = pointed_arch_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], hc=shape.datashape['hc'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'general':
        XY = form.vertices_attributes('xy')
        zub = get_shape_ub_pattern(shape, XY)
        zlb = get_shape_lb_pattern(shape, XY)
    else:
        raise Exception

    i = 0

    for key in form.vertices():
        ub_ = float(zub[i])
        lb_ = float(zlb[i])
        x, y, _ = form.vertex_coordinates(key)
        form.vertex_attribute(key, 'ub', value=ub_)
        form.vertex_attribute(key, 'lb', value=lb_)
        i += 1

    return


def apply_envelope_on_xy(form, c=0.5):
    """ Apply an envelope to the FormDiagram in the plan (x, y) by a given distance.

    Parameters
    ----------
    form : ::FormDiagram::
        The input FormDiagram.
    c : float
        The maximum allowed movement of the nodes in ``x`` or ``y``.

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes
    """

    for key, vertex in form.vertex.items():
        form.vertex_attribute(key, 'xmin', vertex.get('x') - c)
        form.vertex_attribute(key, 'xmax', vertex.get('x') + c)
        form.vertex_attribute(key, 'ymin', vertex.get('y') - c)
        form.vertex_attribute(key, 'ymax', vertex.get('y') + c)

    return


def apply_bounds_on_q(form, qmin=-1e+4, qmax=1e-8):  # Convention compression negative
    """ Apply bounds on the magnitude of the edges'force densities.

    Parameters
    ----------
    form : ::FormDiagram::
        The input FormDiagram
    qmin : float
        The minimum allowed force density ``qmin``.
    qmax : float
        The maximum allowed force density ``qmax``.

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes.
    """

    if isinstance(qmin, list):
        for i, (u, v) in enumerate(form.edges_where({'_is_edge': True})):
            form.edge_attribute((u, v), 'qmin', qmin[i])
            form.edge_attribute((u, v), 'qmax', qmax[i])
    else:
        for u, v in form.edges_where({'_is_edge': True}):
            form.edge_attribute((u, v), 'qmin', qmin)
            form.edge_attribute((u, v), 'qmax', qmax)

    return


def project_mesh_to_middle(mesh, shape=None):
    """ Project a mesh to the middle surface of the shape.

    Parameters
    ----------
    meesh : ::Mesh::
        The input Mesh.
    shape : ::Shape::
        The input Shape with a middle surface.

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes.
    """

    x = mesh.vertices_attribute('x')
    y = mesh.vertices_attribute('y')

    if shape.datashape['type'] == 'dome':
        zt = dome_zt_update(x, y, shape.datashape['radius'], shape.datashape['t'], shape.datashape['center'])
    elif shape.datashape['type'] == 'crossvault':
        zt = crossvault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        zt = pointed_vault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'],
                                         hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
    else:
        XY = mesh.vertices_attributes('xy')
        zt = get_shape_middle_pattern(shape, XY)

    i = 0

    for key in mesh.vertices():
        z = float(zt[i])
        x, y, _ = mesh.vertex_coordinates(key)
        mesh.vertex_attribute(key, 'z', value=z)
        i += 1

    return
