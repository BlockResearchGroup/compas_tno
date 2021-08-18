from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_ub_lb_update
from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.pointed_arch import pointed_arch_ub_lb_update

from compas_tno.utilities.interpolation import get_shape_ub_pattern
from compas_tno.utilities.interpolation import get_shape_lb_pattern


__all__ = [
    'apply_envelope_from_shape',
    'apply_envelope_on_xy',
    'apply_bounds_on_q',

    'apply_envelope_from_shape_proxy',
    'apply_bounds_on_q_proxy',
]


def apply_envelope_from_shape_proxy(formdata, shapedata):  # this works only forlibrary shapes
    # TODO: crate a proper to_data from_data for shapes, and make it happen.

    from compas_tno.diagrams import FormDiagram
    from compas_tno.shapes import Shape

    form = FormDiagram.from_data(formdata)
    shape = Shape.from_library(shapedata)

    apply_envelope_from_shape(form, shape)

    return form.to_data()


def apply_bounds_on_q_proxy(formdata, qmin=-1e+4, qmax=1e-8):  # no need of proxy - change in future

    from compas_tno.diagrams import FormDiagram

    form = FormDiagram.from_data(formdata)

    apply_bounds_on_q(form, qmin=qmin, qmax=qmax)

    return form.to_data()


def apply_envelope_from_shape(form, shape):

    x = form.vertices_attribute('x')  # check if array is necessary here
    y = form.vertices_attribute('y')

    if shape.datashape['type'] == 'dome':
        zub, zlb = dome_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['center'], shape.datashape['radius'])
    elif shape.datashape['type'] == 'crossvault':
        zub, zlb = crossvault_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        zub, zlb = pointed_vault_ub_lb_update(x, y, shape.datashape['thk'], shape.datashape['t'], shape.datashape['xy_span'], hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
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

    for key, vertex in form.vertex.items():
        form.vertex_attribute(key, 'xmin', vertex.get('x') - c)
        form.vertex_attribute(key, 'xmax', vertex.get('x') + c)
        form.vertex_attribute(key, 'ymin', vertex.get('y') - c)
        form.vertex_attribute(key, 'ymax', vertex.get('y') + c)

    return


def apply_bounds_on_q(form, qmin=-1e+4, qmax=1e-8):  # Convention compression negative

    if isinstance(qmin, list):
        for i, (u, v) in enumerate(form.edges_where({'_is_edge': True})):
            form.edge_attribute((u, v), 'qmin', qmin[i])
            form.edge_attribute((u, v), 'qmax', qmax[i])
    else:
        for u, v in form.edges_where({'_is_edge': True}):
            form.edge_attribute((u, v), 'qmin', qmin)
            form.edge_attribute((u, v), 'qmax', qmax)

    return
