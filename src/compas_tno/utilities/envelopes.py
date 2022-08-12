import math

from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.dome import dome_zt_update

from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_middle_update

from compas_tno.shapes.pointed_crossvault import pointed_vault_ub_lb_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update

from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.pointed_arch import pointed_arch_ub_lb_update

from compas_tno.shapes.pavillionvault import pavillionvault_ub_lb_update
from compas_tno.shapes.pavillionvault import pavillionvault_middle_update

from compas_tno.utilities.interpolation import get_shape_ub_pattern
from compas_tno.utilities.interpolation import get_shape_lb_pattern
from compas_tno.utilities.interpolation import get_shape_middle_pattern


__all__ = [
    'apply_envelope_from_shape',
    'apply_envelope_on_xy',
    'apply_envelope_on_xy_from_base',
    'apply_bounds_on_q',
    'project_mesh_to_middle',
    'modify_shapedata_with_spr_angle',
    'apply_envelope_from_shape_proxy',
    'apply_bounds_tub_tlb',
    'apply_bounds_reactions'
]


def apply_envelope_from_shape_proxy(formdata, shapedata):
    """ Apply an envelope (intrados and extrados) to the FormDiagram based on the input shape data

    Parameters
    ----------
    formdata : dict
        The data of the form diagram
    shapedata : dict
        The data of the shape

    Returns
    ----------
    formdata : dict
        The data of the form diagram after assigning the envelope
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
    form : FormDiagram
        The input FormDiagram.
    shape : Shape
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
    elif shape.datashape['type'] == 'pavillionvault':
        zub, zlb = pavillionvault_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['xy_span'], shape.datashape['spr_angle'])
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
    form : FormDiagram
        The input FormDiagram.
    c : float, optional
        The maximum allowed movement of the nodes in ``x`` or ``y``, by default 0.5

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


def apply_envelope_on_xy_from_base(form, form_base, c=0.5):
    """ Apply an envelope to the FormDiagram considering a given distance applied to a base form diagram.

    Parameters
    ----------
    form : FormDiagram
        The input FormDiagram.
    form_base ::FormDiagram
        The base FormDiagram to consider.
    c : float, optional
        The maximum allowed movement of the nodes in ``x`` or ``y``, by default 0.5

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes
    """

    for key, vertex in form_base.vertex.items():
        form.vertex_attribute(key, 'xmin', vertex.get('x') - c)
        form.vertex_attribute(key, 'xmax', vertex.get('x') + c)
        form.vertex_attribute(key, 'ymin', vertex.get('y') - c)
        form.vertex_attribute(key, 'ymax', vertex.get('y') + c)

    return


def apply_bounds_on_q(form, qmin=-1e+4, qmax=1e-8):
    """ Apply bounds on the magnitude of the edges'force densities.

    Parameters
    ----------
    form : FormDiagram
        The input FormDiagram
    qmin : float, optional
        The minimum allowed force density ``qmin``, by default -1e+4
    qmax : float, optional
        The maximum allowed force density ``qmax``, by default 1e-8

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


def apply_bounds_reactions(form, shape, assume_shape=None):
    """ Apply bounds on the magnitude of the allowed increase in thickness of the upper-bound (tub), lower-bound (tlb), and of the reaction vector (tub_reacmax).

    Parameters
    ----------
    form : FormDiagram
        The input FormDiagram
    shape : Shape
        The shape of masonry to constraint the form diagram to.
    assume_shape : dict, optional
        Whether or not consider the settings of a different shape on the process, by default None

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes.
    """

    if assume_shape:
        data = assume_shape
    else:
        data = shape.datashape

    thk = data['thk']

    if data['type'] == 'dome' or data['type'] == 'dome_polar':
        [x0, y0] = data['center'][:2]
        for key in form.vertices_where({'is_fixed': True}):
            x, y, _ = form.vertex_coordinates(key)
            theta = math.atan2((y - y0), (x - x0))
            x_ = thk/2*math.cos(theta)
            y_ = thk/2*math.sin(theta)
            form.vertex_attribute(key, 'b', [x_, y_])

    b_manual = data.get('b_manual', None)
    if data['type'] == 'arch':
        H = data['H']
        L = data['L']
        thk = data['thk']
        radius = H / 2 + (L**2 / (8 * H))
        zc = radius - H
        re = radius + thk/2
        x = math.sqrt(re**2 - zc**2)
        for key in form.vertices_where({'is_fixed': True}):
            form.vertex_attribute(key, 'b', [x - L/2, 0.0])
            if b_manual:
                form.vertex_attribute(key, 'b', [b_manual, 0.0])
                print('Applied b manual')

    if data['type'] == 'dome_spr':
        x0 = data['center'][0]
        y0 = data['center'][1]
        radius = data['radius']
        [_, theta_f] = data['theta']
        r_proj_e = (radius + thk/2) * math.sin(theta_f)
        r_proj_m = (radius) * math.sin(theta_f)
        delt = r_proj_e - r_proj_m
        for key in form.vertices_where({'is_fixed': True}):
            x, y, _ = form.vertex_coordinates(key)
            theta = math.atan2((y - y0), (x - x0))
            x_ = delt*math.cos(theta)
            y_ = delt*math.sin(theta)
            form.vertex_attribute(key, 'b', [x_, y_])

    if data['type'] == 'pavillionvault':
        x0, x1 = data['xy_span'][0]
        y0, y1 = data['xy_span'][1]
        for key in form.vertices_where({'is_fixed': True}):
            x, y, _ = form.vertex_coordinates(key)
            if x == x0:
                form.vertex_attribute(key, 'b', [-thk/2, 0])
            elif x == x1:
                form.vertex_attribute(key, 'b', [+thk/2, 0])
            if y == y0:
                if form.vertex_attribute(key, 'b'):
                    b = form.vertex_attribute(key, 'b')
                    b[1] = -thk/2
                    form.vertex_attribute(key, 'b', b)
                else:
                    form.vertex_attribute(key, 'b', [0, -thk/2])
            elif y == y1:
                if form.vertex_attribute(key, 'b'):
                    b = form.vertex_attribute(key, 'b')
                    b[1] = +thk/2
                    form.vertex_attribute(key, 'b', b)
                else:
                    form.vertex_attribute(key, 'b', [0, +thk/2])

    return


def apply_bounds_tub_tlb(form, tubmax=0.5, tlbmax=0.5):
    """ Apply bounds on the magnitude of the allowed increase in thickness of the upper-bound (tub), lower-bound (tlb), and of the reaction vector (tub_reacmax).

    Parameters
    ----------
    form : FormDiagram
        The input FormDiagram
    tubmax : float, optional
        The maximum increase in thickness of the extrados. The default value is ``0.5``.
    tlbmax : float, optional
        The maximum increase in thickness of the intrados. The default value is ``0.5``.

    Returns
    ----------
    None
        The formdiagram is updated in place in the attributes.
    """

    for vertex in form.vertices():
        form.vertex_attribute(vertex, 'tubmax', tubmax)
        form.vertex_attribute(vertex, 'tlbmax', tlbmax)

    return


def project_mesh_to_middle(mesh, shape=None):
    """ Project a mesh to the middle surface of the shape.

    Parameters
    ----------
    mesh : Mesh
        The input Mesh.
    shape : Shape
        The input Shape with a middle surface, by default None.

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
    elif shape.datashape['type'] == 'pavillionvault':
        zt = pavillionvault_middle_update(x, y,  shape.datashape['t'],  xy_span=shape.datashape['xy_span'], spr_angle=shape.datashape['spr_angle'])
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


def modify_shapedata_with_spr_angle(datashape):
    """ Modify the Shape data to account for an springing angle.
    This is done by increasing the ``xy_span`` limits and works exclusivvely for rectandular diagrams.

    Parameters
    ----------
    datashape : dict
        The data to construct the Shape.

    Returns
    ----------
    datashape : dict
        The updated datashape.
    """

    typevault = datashape.get('type')
    spr_angle = datashape.get('spr_angle', 0.0)

    if spr_angle:
        if typevault in ['crossvault', 'pointed_crossvault', 'pavillionvault', ':parabolic_shell', 'domicalvault']:
            xy_span = datashape['xy_span']
            span_x = xy_span[0][1] - xy_span[0][0]
            span_y = xy_span[1][1] - xy_span[1][0]

            A = 1/math.cos(math.radians(spr_angle))
            print('deg, discretisation:', A, spr_angle)
            xy_span_enlarged = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]

            datashape['xy_span'] = xy_span_enlarged
        else:
            raise NotImplementedError

    return datashape
