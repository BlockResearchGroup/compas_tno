from typing import List

import numpy.typing as npt

from compas_tna.envelope import Envelope

# TODO: Add the case for the general MeshEnvelope
# from compas_tno.shapes import general_b_update_with_n
# from compas_tno.shapes import general_db_with_n
# from compas_tno.shapes import general_db_with_t_intrados
# from compas_tno.shapes import general_db_with_t_middle_constant
# from compas_tno.shapes import general_db_with_t_middle_variable
# from compas_tno.shapes import general_dub_dlb_with_n
# from compas_tno.shapes import general_dub_dlb_with_s
# from compas_tno.shapes import general_ub_lb_update_with_n
# from compas_tno.shapes import general_ub_lb_update_with_s
# from compas_tno.shapes import general_ub_lb_update_with_t_intrados
# from compas_tno.shapes import general_ub_lb_update_with_t_middle_constant
# from compas_tno.shapes import general_ub_lb_update_with_t_middle_variable


def ub_lb_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    t: float,
    envelope: Envelope,
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
    envelope : Envelope
        Envelope object to be updated
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

    if envelope.is_parametric:
        return envelope.compute_bounds(x, y, thk)
    else:
        raise NotImplementedError("Future implementation")

    #     if "t" in variables:
    #         thickness_type = envelope.thickness_type

    #         if thickness_type == "constant":
    #             return general_ub_lb_update_with_t_middle_constant(thk, s, shape.middle, t)  # here open up something
    #         if thickness_type == "variable":
    #             return general_ub_lb_update_with_t_middle_variable(thk, s, shape.middle, t)  # here open up something
    #         if thickness_type == "intrados":
    #             return general_ub_lb_update_with_t_intrados(thk, lb, shape.intrados, t)  # here open up something

    #     elif "s" in variables:
    #         return general_ub_lb_update_with_s(ub, lb, thk)  # thk is 's' in this equation

    #     elif "n" in variables:
    #         return general_ub_lb_update_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation

    # raise Exception


def dub_dlb_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    t: float,
    envelope: Envelope,
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
    envelope : Envelope
        Envelope object to be updated
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

    if envelope.is_parametric:
        return envelope.compute_bounds_derivatives(x, y, thk)
    else:
        raise NotImplementedError("Future implementation")

        # TODO: Add the case for the general MeshEnvelope

        # if "t" in variables:
        #     thickness_type = shape.parameters["thickness_type"]

        #     if thickness_type == "constant":
        #         return general_db_with_t_middle_constant(s, shape.middle)  # here open up something
        #     if thickness_type == "variable":
        #         return general_db_with_t_middle_variable(s, shape.middle)  # here open up something
        #     if thickness_type == "intrados":
        #         return general_db_with_t_intrados(lb, shape.intrados)  # here open up something

        # elif "s" in variables:
        #     return general_dub_dlb_with_s(ub, lb)  # thk is 's' in this equation
        # elif "n" in variables:
        #     return general_dub_dlb_with_n(ub, lb, thk, shape.intrados, shape.extrados, t)  # thk is 'n' in this equation


def b_update(
    x: npt.NDArray,
    y: npt.NDArray,
    thk: float,
    fixed: list,
    envelope: Envelope,
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
    envelope : Envelope
        Envelope object to be updated
    b : array
        Current ``b`` limits
    variables : list
        List with the variables passed

    Returns
    -------
    b : list
        New ``b`` limits
    """

    if envelope.is_parametric:
        return envelope.compute_bound_react(x, y, thk, fixed)
    else:
        raise NotImplementedError("Future implementation")

    # TODO: Add the case for the general MeshEnvelope
    #     if "n" in variables:
    #         if shape.parameters["base_structure"]["type"] == "dome" or shape.parameters["base_structure"]["type"] == "dome_polar":
    #             return dome_b_update_with_n(x, y, thk, fixed, b, center=shape.parameters["base_structure"]["center"])

    # raise Exception


def db_update(x: npt.NDArray, y: npt.NDArray, thk: float, fixed: List[int], envelope: Envelope, b: npt.NDArray, variables: List[str]):
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
    envelope : Envelope
        Envelope object to be updated
    b : array
        Current ``b`` limits
    variables : list
        List with the variables passed

    Returns
    -------
    db
        Sensitivities of the ``b`` limits in the point
    """

    if envelope.is_parametric:
        return envelope.compute_bound_react_derivatives(x, y, thk, fixed)
    else:
        raise NotImplementedError("Future implementation")

    # TODO: Add the case for the general MeshEnvelope
    # if shape.parameters["type"] == "arch":
    #     return arch_db(x, y, thk, fixed, H=shape.parameters["H"], L=shape.parameters["L"], x0=shape.parameters["x0"])

    # if shape.parameters["type"] == "dome":
    #     return dome_db(x, y, thk, fixed, center=shape.parameters["center"], radius=shape.parameters["radius"])

    # if shape.parameters["type"] == "pavillionvault":
    #     return pavillionvault_db(x, y, thk, fixed, xy_span=shape.parameters["xy_span"])

    # if shape.parameters["type"] == "general":
    #     if "n" in variables:
    #         if shape.parameters["base_structure"]["type"] == "dome":
    #             return dome_db_with_n(x, y, fixed, center=shape.parameters["base_structure"]["center"])

    # raise Exception
