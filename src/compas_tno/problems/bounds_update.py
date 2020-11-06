from compas_tno.shapes.circular_arch import arch_ub_lb_update
from compas_tno.shapes.circular_arch import arch_b_update
from compas_tno.shapes.circular_arch import arch_dub_dlb
from compas_tno.shapes.circular_arch import arch_db

from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.dome import dome_b_update
from compas_tno.shapes.dome import dome_dub_dlb
from compas_tno.shapes.dome import dome_db
from compas_tno.shapes.dome import dome_b_update_with_n
from compas_tno.shapes.dome import dome_db_with_n

from compas_tno.shapes.pavillionvault import pavillionvault_ub_lb_update
from compas_tno.shapes.pavillionvault import pavillionvault_b_update
from compas_tno.shapes.pavillionvault import pavillionvault_dub_dlb
from compas_tno.shapes.pavillionvault import pavillionvault_db

from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_dub_dlb

from compas_tno.shapes.general import general_ub_lb_update_with_s
from compas_tno.shapes.general import general_dub_dlb_with_s
from compas_tno.shapes.general import general_ub_lb_update_with_n
from compas_tno.shapes.general import general_dub_dlb_with_n
from compas_tno.shapes.general import general_b_update_with_n
from compas_tno.shapes.general import general_db_with_n


def ub_lb_update(x, y, thk, t, shape, ub, lb, variables):

    if shape.data['type'] == 'arch':
        return arch_ub_lb_update(x, y, thk, t, H=shape.data['H'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'dome':
        return dome_ub_lb_update(x, y, thk, t, center=shape.data['center'], radius=shape.data['radius'])
    elif shape.data['type'] == 'crossvault':
        return crossvault_ub_lb_update(x, y, thk, t, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'pavillionvault':
        return pavillionvault_ub_lb_update(x, y, thk, t, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'general':
        if 's' in variables:
            return general_ub_lb_update_with_s(ub, lb, thk)  # thk is 's' in this equation
        elif 'n' in variables:
            return general_ub_lb_update_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def dub_dlb_update(x, y, thk, t, shape, ub, lb, variables):

    if shape.data['type'] == 'arch':
        return arch_dub_dlb(x, y, thk, t, H=shape.data['H'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'dome':
        return dome_dub_dlb(x, y, thk, t, center=shape.data['center'], radius=shape.data['radius'])
    elif shape.data['type'] == 'crossvault':
        return crossvault_dub_dlb(x, y, thk, t, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'pavillionvault':
        return pavillionvault_dub_dlb(x, y, thk, t, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'general':
        if 's' in variables:
            return general_dub_dlb_with_s(ub, lb)  # thk is 's' in this equation
        elif 'n' in variables:
            return general_dub_dlb_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def b_update(x, y, thk, fixed, shape, b, variables):

    if shape.data['type'] == 'arch':
        return arch_b_update(x, y, thk, fixed, H=shape.data['H'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'dome':
        return dome_b_update(x, y, thk, fixed, center=shape.data['center'], radius=shape.data['radius'])
    elif shape.data['type'] == 'pavillionvault':
        return pavillionvault_b_update(x, y, thk, fixed, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'general':
        if 'n' in variables:
            if shape.data['base_structure']['type'] == 'dome' or shape.data['base_structure']['type'] == 'dome_polar':
                return dome_b_update_with_n(x, y, thk, fixed, b, center=shape.data['base_structure']['center'])
            else:
                return general_b_update_with_n(b, thk, fixed)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def db_update(x, y, thk, fixed, shape, b, variables):

    if shape.data['type'] == 'arch':
        return arch_db(x, y, thk, fixed, H=shape.data['H'], L=shape.data['L'], x0=shape.data['x0'])
    elif shape.data['type'] == 'dome':
        return dome_db(x, y, thk, fixed, center=shape.data['center'], radius=shape.data['radius'])
    elif shape.data['type'] == 'pavillionvault':
        return pavillionvault_db(x, y, thk, fixed, xy_span=shape.data['xy_span'])
    elif shape.data['type'] == 'general':
        if 'n' in variables:
            if shape.data['base_structure']['type'] == 'dome':
                return dome_db_with_n(x, y, fixed, center=shape.data['base_structure']['center'])
            else:
                return general_b_update_with_n(b, thk, fixed)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception
