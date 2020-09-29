from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.circular_arch import arch_b_update
from compas_tno.shapes.circular_arch import arch_dub_dlb
from compas_tno.shapes.circular_arch import arch_db

from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.dome import dome_b_update
from compas_tno.shapes.dome import dome_dub_dlb
from compas_tno.shapes.dome import dome_db

from compas_tno.shapes.pavillionvault import pavillionvault_ub_lb_update
from compas_tno.shapes.pavillionvault import pavillionvault_b_update
from compas_tno.shapes.pavillionvault import pavillionvault_dub_dlb
from compas_tno.shapes.pavillionvault import pavillionvault_db

from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_dub_dlb


def ub_lb_update(x, y, thk, t, shape_data):

    if shape_data['type'] == 'arch':
        return arch_ub_lb_update(x, y, thk, t, H=shape_data['H'], L=shape_data['L'], x0=shape_data['x0'])
    elif shape_data['type'] == 'dome':
        return dome_ub_lb_update(x, y, thk, t, center=shape_data['center'], radius=shape_data['radius'])
    elif shape_data['type'] == 'crossvault':
        return crossvault_ub_lb_update(x, y, thk, t, xy_span=shape_data['xy_span'])
    elif shape_data['type'] == 'pavillionvault':
        return pavillionvault_ub_lb_update(x, y, thk, t, xy_span=shape_data['xy_span'])
    else:
        raise Exception


def dub_dlb_update(x, y, thk, t, shape_data):

    if shape_data['type'] == 'arch':
        return arch_dub_dlb(x, y, thk, t, H=shape_data['H'], L=shape_data['L'], x0=shape_data['x0'])
    elif shape_data['type'] == 'dome':
        return dome_dub_dlb(x, y, thk, t, center=shape_data['center'], radius=shape_data['radius'])
    elif shape_data['type'] == 'crossvault':
        return crossvault_dub_dlb(x, y, thk, t, xy_span=shape_data['xy_span'])
    elif shape_data['type'] == 'pavillionvault':
        return pavillionvault_dub_dlb(x, y, thk, t, xy_span=shape_data['xy_span'])
    else:
        raise Exception


def b_update(x, y, thk, fixed, shape_data):

    if shape_data['type'] == 'arch':
        return arch_b_update(x, y, thk, fixed, H=shape_data['H'], L=shape_data['L'], x0=shape_data['x0'])
    elif shape_data['type'] == 'dome':
        return dome_b_update(x, y, thk, fixed, center=shape_data['center'], radius=shape_data['radius'])
    elif shape_data['type'] == 'pavillionvault':
        return pavillionvault_b_update(x, y, thk, fixed, xy_span=shape_data['xy_span'])
    else:
        raise Exception


def db_update(x, y, thk, fixed, shape_data):

    if shape_data['type'] == 'arch':
        return arch_db(x, y, thk, fixed, H=shape_data['H'], L=shape_data['L'], x0=shape_data['x0'])
    elif shape_data['type'] == 'dome':
        return dome_db(x, y, thk, fixed, center=shape_data['center'], radius=shape_data['radius'])
    elif shape_data['type'] == 'pavillionvault':
        return pavillionvault_db(x, y, thk, fixed, xy_span=shape_data['xy_span'])
    else:
        raise Exception
