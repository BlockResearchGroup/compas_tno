import numpy.typing as npt
from typing import List, Tuple,TYPE_CHECKING
if TYPE_CHECKING:
    from compas_tno.problems import Problem

# TODO: This needs to be taken care by the SurfaceModel
# Since the Model will be created from a template, the bounds will be updated based on the template
# and the thickness will be updated based on the template

from compas_tno.shapes import Shape
from compas_tno.shapes import arch_b_update
from compas_tno.shapes import arch_db
from compas_tno.shapes import arch_dub_dlb
from compas_tno.shapes import arch_ub_lb_update
from compas_tno.shapes import crossvault_dub_dlb
from compas_tno.shapes import crossvault_ub_lb_update
from compas_tno.shapes import dome_b_update
from compas_tno.shapes import dome_b_update_with_n
from compas_tno.shapes import dome_db
from compas_tno.shapes import dome_db_with_n
from compas_tno.shapes import dome_dub_dlb
from compas_tno.shapes import dome_ub_lb_update
from compas_tno.shapes import general_db_with_t_intrados
from compas_tno.shapes import general_db_with_t_middle_constant
from compas_tno.shapes import general_db_with_t_middle_variable
from compas_tno.shapes import general_dub_dlb_with_n
from compas_tno.shapes import general_dub_dlb_with_s
from compas_tno.shapes import general_ub_lb_update_with_n
from compas_tno.shapes import general_ub_lb_update_with_s
from compas_tno.shapes import general_ub_lb_update_with_t_intrados

# from compas_tno.shapes import general_b_update_with_n
# from compas_tno.shapes import general_db_with_n
from compas_tno.shapes import general_ub_lb_update_with_t_middle_constant
from compas_tno.shapes import general_ub_lb_update_with_t_middle_variable
from compas_tno.shapes import pavillionvault_b_update
from compas_tno.shapes import pavillionvault_db
from compas_tno.shapes import pavillionvault_dub_dlb
from compas_tno.shapes import pavillionvault_ub_lb_update
from compas_tno.shapes import pointed_arch_dub_dlb
from compas_tno.shapes import pointed_arch_ub_lb_update
from compas_tno.shapes import pointed_vault_dub_dlb
from compas_tno.shapes import pointed_vault_ub_lb_update


def ub_lb_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    t: float,
    shape: Shape,
    ub: npt.NDArray,
    lb: npt.NDArray,
    s: npt.NDArray,
    variables: list,
) -> tuple[npt.NDArray, npt.NDArray]:
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

    if shape.parameters["type"] == "arch":
        return arch_ub_lb_update(x, y, thk, t, H=shape.parameters["H"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "pointed_arch":
        return pointed_arch_ub_lb_update(x, y, thk, t, hc=shape.parameters["hc"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "dome":
        return dome_ub_lb_update(x, y, thk, t, center=shape.parameters["center"], radius=shape.parameters["radius"])

    if shape.parameters["type"] == "crossvault":
        return crossvault_ub_lb_update(x, y, thk, t, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "pavillionvault":
        return pavillionvault_ub_lb_update(x, y, thk, t, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "pointed_crossvault":
        return pointed_vault_ub_lb_update(x, y, thk, t, xy_span=shape.parameters["xy_span"], hc=shape.parameters["hc"], he=shape.parameters["he"], hm=shape.parameters["hm"])

    if shape.parameters["type"] == "general":
        if "t" in variables:
            thickness_type = shape.parameters["thickness_type"]

            if thickness_type == "constant":
                return general_ub_lb_update_with_t_middle_constant(thk, s, shape.middle, t)  # here open up something
            if thickness_type == "variable":
                return general_ub_lb_update_with_t_middle_variable(thk, s, shape.middle, t)  # here open up something
            if thickness_type == "intrados":
                return general_ub_lb_update_with_t_intrados(thk, lb, shape.intrados, t)  # here open up something

        elif "s" in variables:
            return general_ub_lb_update_with_s(ub, lb, thk)  # thk is 's' in this equation

        elif "n" in variables:
            return general_ub_lb_update_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation

    raise Exception


def dub_dlb_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    t: float,
    shape: Shape,
    ub: npt.NDArray,
    lb: npt.NDArray,
    s: npt.NDArray,
    variables: list,
) -> tuple[npt.NDArray, npt.NDArray]:
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
    s : array
        Original level of the supports
    variables : list
        List with the variables passed

    Returns
    -------
    dub, dlb
        New sensitivities of the bounds in the point analysed

    """

    if shape.parameters["type"] == "arch":
        return arch_dub_dlb(x, y, thk, t, H=shape.parameters["H"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "pointed_arch":
        return pointed_arch_dub_dlb(x, y, thk, t, hc=shape.parameters["hc"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "dome":
        return dome_dub_dlb(x, y, thk, t, center=shape.parameters["center"], radius=shape.parameters["radius"])

    if shape.parameters["type"] == "crossvault":
        return crossvault_dub_dlb(x, y, thk, t, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "pavillionvault":
        return pavillionvault_dub_dlb(x, y, thk, t, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "pointed_crossvault":
        return pointed_vault_dub_dlb(x, y, thk, t, xy_span=shape.parameters["xy_span"], hc=shape.parameters["hc"], he=shape.parameters["he"], hm=shape.parameters["hm"])

    if shape.parameters["type"] == "general":
        if "t" in variables:
            thickness_type = shape.parameters["thickness_type"]

            if thickness_type == "constant":
                return general_db_with_t_middle_constant(s, shape.middle)  # here open up something
            if thickness_type == "variable":
                return general_db_with_t_middle_variable(s, shape.middle)  # here open up something
            if thickness_type == "intrados":
                return general_db_with_t_intrados(lb, shape.intrados)  # here open up something

        elif "s" in variables:
            return general_dub_dlb_with_s(ub, lb)  # thk is 's' in this equation
        elif "n" in variables:
            return general_dub_dlb_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation

    raise Exception


def b_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    fixed: list,
    shape: Shape,
    b: npt.NDArray,
    variables: list,
) -> list:
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
    variables : list
        List with the variables passed

    Returns
    -------
    b : list
        New ``b`` limits
    """

    if shape.parameters["type"] == "arch":
        return arch_b_update(x, y, thk, fixed, H=shape.parameters["H"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "dome":
        return dome_b_update(x, y, thk, fixed, center=shape.parameters["center"], radius=shape.parameters["radius"])

    if shape.parameters["type"] == "pavillionvault":
        return pavillionvault_b_update(x, y, thk, fixed, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "general":
        if "n" in variables:
            if shape.parameters["base_structure"]["type"] == "dome" or shape.parameters["base_structure"]["type"] == "dome_polar":
                return dome_b_update_with_n(x, y, thk, fixed, b, center=shape.parameters["base_structure"]["center"])

    raise Exception


def db_update(x: npt.NDArray, y: npt.NDArray, thk: float, fixed: List[int], shape: Shape, b: npt.NDArray, variables: List[str]):
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
    variables : list
        List with the variables passed

    Returns
    -------
    db
        Sensitivities of the ``b`` limits in the point
    """

    if shape.parameters["type"] == "arch":
        return arch_db(x, y, thk, fixed, H=shape.parameters["H"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    if shape.parameters["type"] == "dome":
        return dome_db(x, y, thk, fixed, center=shape.parameters["center"], radius=shape.parameters["radius"])

    if shape.parameters["type"] == "pavillionvault":
        return pavillionvault_db(x, y, thk, fixed, xy_span=shape.parameters["xy_span"])

    if shape.parameters["type"] == "general":
        if "n" in variables:
            if shape.parameters["base_structure"]["type"] == "dome":
                return dome_db_with_n(x, y, fixed, center=shape.parameters["base_structure"]["center"])

    raise Exception
