from compas_tno.shapes import arch_ub_lb_update
from compas_tno.shapes import arch_b_update
from compas_tno.shapes import arch_dub_dlb
from compas_tno.shapes import arch_db

from compas_tno.shapes import pointed_arch_ub_lb_update
from compas_tno.shapes import pointed_arch_dub_dlb

from compas_tno.shapes import dome_ub_lb_update
from compas_tno.shapes import dome_b_update
from compas_tno.shapes import dome_dub_dlb
from compas_tno.shapes import dome_db
from compas_tno.shapes import dome_b_update_with_n
from compas_tno.shapes import dome_db_with_n

from compas_tno.shapes import pavillionvault_ub_lb_update
from compas_tno.shapes import pavillionvault_b_update
from compas_tno.shapes import pavillionvault_dub_dlb
from compas_tno.shapes import pavillionvault_db

from compas_tno.shapes import crossvault_ub_lb_update
from compas_tno.shapes import crossvault_dub_dlb

from compas_tno.shapes import pointed_vault_ub_lb_update
from compas_tno.shapes import pointed_vault_dub_dlb

from compas_tno.shapes import general_ub_lb_update_with_s
from compas_tno.shapes import general_dub_dlb_with_s
from compas_tno.shapes import general_ub_lb_update_with_n
from compas_tno.shapes import general_dub_dlb_with_n
# from compas_tno.shapes import general_b_update_with_n
# from compas_tno.shapes import general_db_with_n
from compas_tno.shapes import general_ub_lb_update_with_t_middle_constant
from compas_tno.shapes import general_db_with_t_middle_constant
from compas_tno.shapes import general_ub_lb_update_with_t_middle_variable
from compas_tno.shapes import general_db_with_t_middle_variable
from compas_tno.shapes import general_ub_lb_update_with_t_intrados
from compas_tno.shapes import general_db_with_t_intrados


def ub_lb_update(x, y, thk, t, shape, ub, lb, s, variables):
    """Function to update the ub-lb vertical bounds of the vertices.

    Parameters
    ----------
    x : array
        xcoordinates of the pattern
    y : array
        y-coordinates of the pattern
    thk : float
        thickness to update the bounds to
    t : float
        thickness considered to the outer vertices
    shape : Shape
        Shape object to be updated
    ub : array
        Current upper-bound limits
    lb : array
        Current lower-bound limits
    s : aray
        Original level of the supports
    variables : list
        List with the variables passed

    Returns
    -------
    ub, lb: array
        New position of the bounds in the point analysed
    """

    if shape.datashape['type'] == 'arch':
        return arch_ub_lb_update(x, y, thk, t, H=shape.datashape['H'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    if shape.datashape['type'] == 'pointed_arch':
        return pointed_arch_ub_lb_update(x, y, thk, t, hc=shape.datashape['hc'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'dome':
        return dome_ub_lb_update(x, y, thk, t, center=shape.datashape['center'], radius=shape.datashape['radius'])
    elif shape.datashape['type'] == 'crossvault':
        return crossvault_ub_lb_update(x, y, thk, t, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pavillionvault':
        return pavillionvault_ub_lb_update(x, y, thk, t, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        return pointed_vault_ub_lb_update(x, y, thk, t, xy_span=shape.datashape['xy_span'], hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
    elif shape.datashape['type'] == 'general':
        if 't' in variables:
            thickness_type = shape.datashape['thickness_type']
            if thickness_type == 'constant':
                return general_ub_lb_update_with_t_middle_constant(thk, s, shape.middle, t)  # here open up something
            elif thickness_type == 'variable':
                return general_ub_lb_update_with_t_middle_variable(thk, s, shape.middle, t)  # here open up something
            elif thickness_type == 'intrados':
                return general_ub_lb_update_with_t_intrados(thk, lb, shape.intrados, t)  # here open up something
            else:
                raise Exception
        elif 's' in variables:
            return general_ub_lb_update_with_s(ub, lb, thk)  # thk is 's' in this equation
        elif 'n' in variables:
            return general_ub_lb_update_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def dub_dlb_update(x, y, thk, t, shape, ub, lb, s, variables):
    """Function to update the derivatives of the ub-lb vertical bounds of the vertices.

    Parameters
    ----------
    x : array
        xcoordinates of the pattern
    y : array
        y-coordinates of the pattern
    thk : float
        thickness to update the bounds to
    t : float
        thickness considered to the outer vertices
    shape : Shape
        Shape object to be updated
    ub : array
        Current upper-bound limits
    lb : array
        Current lower-bound limits
    s : aray
        Original level of the supports
    variables : list
        List with the variables passed

    Returns
    -------
    dub, dlb
        New sensitivities of the bounds in the point analysed
    """

    if shape.datashape['type'] == 'arch':
        return arch_dub_dlb(x, y, thk, t, H=shape.datashape['H'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    if shape.datashape['type'] == 'pointed_arch':
        return pointed_arch_dub_dlb(x, y, thk, t, hc=shape.datashape['hc'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'dome':
        return dome_dub_dlb(x, y, thk, t, center=shape.datashape['center'], radius=shape.datashape['radius'])
    elif shape.datashape['type'] == 'crossvault':
        return crossvault_dub_dlb(x, y, thk, t, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pavillionvault':
        return pavillionvault_dub_dlb(x, y, thk, t, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'pointed_crossvault':
        return pointed_vault_dub_dlb(x, y, thk, t, xy_span=shape.datashape['xy_span'], hc=shape.datashape['hc'], he=shape.datashape['he'], hm=shape.datashape['hm'])
    elif shape.datashape['type'] == 'general':
        if 't' in variables:
            thickness_type = shape.datashape['thickness_type']
            if thickness_type == 'constant':
                return general_db_with_t_middle_constant(s, shape.middle)  # here open up something
            elif thickness_type == 'variable':
                return general_db_with_t_middle_variable(s, shape.middle)  # here open up something
            elif thickness_type == 'intrados':
                return general_db_with_t_intrados(lb, shape.intrados)  # here open up something
            else:
                raise Exception
        elif 's' in variables:
            return general_dub_dlb_with_s(ub, lb)  # thk is 's' in this equation
        elif 'n' in variables:
            return general_dub_dlb_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def b_update(x, y, thk, fixed, shape, b, variables):
    """Function to update the limits of the extension of the reaction forces on the support vertices.

    Parameters
    ----------
    x : array
        xcoordinates of the pattern
    y : array
        y-coordinates of the pattern
    thk : float
        thickness to update the bounds to
    fixed : list
        LIst with the indices of the fixed vertices
    shape : Shape
        Shape object to be updated
    b : array
        Current ``b`` limits

    Returns
    -------
    b : list
        New ``b`` limits
    """

    if shape.datashape['type'] == 'arch':
        return arch_b_update(x, y, thk, fixed, H=shape.datashape['H'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'dome':
        return dome_b_update(x, y, thk, fixed, center=shape.datashape['center'], radius=shape.datashape['radius'])
    elif shape.datashape['type'] == 'pavillionvault':
        return pavillionvault_b_update(x, y, thk, fixed, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'general':
        if 'n' in variables:
            if shape.datashape['base_structure']['type'] == 'dome' or shape.datashape['base_structure']['type'] == 'dome_polar':
                return dome_b_update_with_n(x, y, thk, fixed, b, center=shape.datashape['base_structure']['center'])
            else:
                raise Exception
                # return general_b_update_with_n(b, thk, fixed)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception


def db_update(x, y, thk, fixed, shape, b, variables):
    """Function to update the derrivatives of the limits of the extension of the reaction forces on the support vertices.

    Parameters
    ----------
    x : array
        xcoordinates of the pattern
    y : array
        y-coordinates of the pattern
    thk : float
        thickness to update the bounds to
    fixed : list
        LIst with the indices of the fixed vertices
    shape : Shape
        Shape object to be updated
    b : array
        Current ``b`` limits

    Returns
    -------
    db
        Sensitivities of the ``b`` limits in the point
    """

    if shape.datashape['type'] == 'arch':
        return arch_db(x, y, thk, fixed, H=shape.datashape['H'], L=shape.datashape['L'], x0=shape.datashape['x0'])
    elif shape.datashape['type'] == 'dome':
        return dome_db(x, y, thk, fixed, center=shape.datashape['center'], radius=shape.datashape['radius'])
    elif shape.datashape['type'] == 'pavillionvault':
        return pavillionvault_db(x, y, thk, fixed, xy_span=shape.datashape['xy_span'])
    elif shape.datashape['type'] == 'general':
        if 'n' in variables:
            if shape.datashape['base_structure']['type'] == 'dome':
                return dome_db_with_n(x, y, fixed, center=shape.datashape['base_structure']['center'])
            else:
                raise Exception
                # return general_b_update_with_n(b, thk, fixed)  # thk is 'n' in this equation
        else:
            raise Exception
    else:
        raise Exception
