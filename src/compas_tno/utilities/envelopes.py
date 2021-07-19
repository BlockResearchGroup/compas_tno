from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_ub_lb_update
from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.pointed_arch import pointed_arch_ub_lb_update


__all__ = [
    'apply_envelope_from_shape',
    'apply_envelope_on_xy',
    'apply_bounds_on_q',
]


def apply_envelope_from_shape(form, shape):

    x = form.vertices_attribute('x')  # check if array is necessary here
    y = form.vertices_attribute('y')

    if shape.data['type'] == 'dome':
        zub, zlb = dome_ub_lb_update(x, y, shape.data['thk'], shape.data['t'], shape.data['center'], shape.data['radius'])
    elif shape.data['type'] == 'crossvault':
        zub, zlb = crossvault_ub_lb_update(x, y, shape.data['thk'], shape.data['t'], shape.data['xy_span'])
    elif shape.data['type'] == 'pointed_crossvault':
        zub, zlb = pointed_vault_ub_lb_update(x, y, shape.data['thk'], shape.data['t'], shape.data['xy_span'], hc=shape.data['hc'], he=shape.data['he'], hm=shape.data['hm'])
    elif shape.data['type'] == 'arch':
        zub, zlb = arch_ub_lb_update(x, y, shape.data['thk'], shape.data['t'], H=shape.data['H'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'pointed_arch':
        zub, zlb = pointed_arch_ub_lb_update(x, y, shape.data['thk'], shape.data['t'], hc=shape.data['hc'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'general':
        XY = form.vertices_attributes('xy')
        zub = shape.get_ub_pattern(XY)
        zlb = shape.get_lb_pattern(XY)
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
