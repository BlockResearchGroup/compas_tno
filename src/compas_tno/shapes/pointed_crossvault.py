from numpy import array
from numpy import linspace
from numpy import ones
from numpy import zeros

from compas_tno.shapes import MeshDos
from compas_tno.shapes import rectangular_topology

import math


def pointed_vault_heightfields_proxy(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, discretisation=[20, 20], hc=8.0, he=None, hm=None, tol=10e-6, t=0.0, *args, **kwargs):
    """Set pointed cross vault heights through a proxy

    Parameters
    ----------
    xy_span : list, optional
        [description], by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        [description], by default 0.5
    discretisation : list, optional
        [description], by default [10, 10]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    he : [float, float, float, float], optional
        Height of the opening mid-span for each of the quadrants, by default None
    hm : [float, float, float, float], optional
        Height of each quadrant center (spadrel), by default None
    tol : float, optional
        Tolerance, by default 10e-6
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0

    Returns
    -------
    intradosdata
        The data of MeshDos for the intrados of the shape
    extradosdata
        The data of MeshDos for the extrados of the shape
    middledata
        The data of MeshDos for the middle of the shape

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4
    """
    intrados, extrados, middle = pointed_vault_heightfields(xy_span=xy_span, thk=thk, discretisation=discretisation, hc=hc, he=he, hm=hm, tol=tol, t=t)
    return intrados.to_data(), extrados.to_data(), middle.to_data()


def pointed_vault_heightfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, discretisation=[10, 10], hc=8.0, he=None, hm=None,  t=0.0, tol=0.00):
    """Set pointed cross vault heights

    Parameters
    ----------
    xy_span : list, optional
        [description], by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        [description], by default 0.5
    discretisation : list, optional
        [description], by default [10, 10]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    he : [float, float, float, float], optional
        Height of the opening mid-span for each of the quadrants, by default None
    hm : [float, float, float, float], optional
        Height of each quadrant center (spadrel), by default None
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4
    """

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    density_x = discretisation[0]
    density_y = discretisation[1]
    x = linspace(x0, x1, num=density_x+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0, y1, num=density_y+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)

    xi, yi, faces_i = rectangular_topology(x, y)

    zt = pointed_vault_middle_update(xi, yi, t, xy_span=xy_span, hc=hc, he=he, hm=hm, tol=tol)
    xyzt = array([xi, yi, zt.flatten()]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    zub, zlb = pointed_vault_ub_lb_update(xi, yi, thk, t, xy_span=xy_span, hc=hc, he=he, hm=hm, tol=tol)
    xyzub = array([xi, yi, zub.flatten()]).transpose()
    xyzlb = array([xi, yi, zlb.flatten()]).transpose()
    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)

    return intrados, extrados, middle


def pointed_vault_middle_update(x, y, thk, xy_span=[[0.0, 10.0], [0.0, 10.0]], hc=8.0, he=None, hm=None, tol=1e-6):
    """Update middle of a pointed vault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the arch
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    he : [float, float, float, float], optional
        Height of the opening mid-span for each of the quadrants, by default None
    hm : [float, float, float, float], optional
        Height of each quadrant center (spadrel), by default None
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    middle : array
        Values of the middle surface in the points
    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]
    lx = x1 - x0
    ly = y1 - y0

    if he and hm is None:
        h1, k1, r1 = _circle_3points_xy([x0, he[1]], [(x1+x0)/2, hc], [x1, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = _circle_3points_xy([y0, he[3]], [(y1+y0)/2, hc], [y1, he[2]])
        h4, k4, r4 = h3, k3, r3
    elif hm and he:
        h1, k1, r1 = _circle_3points_xy([(x1+x0)/2, hc], [3*(x1+x0)/4, hm[0]], [x1, he[0]])
        h2, k2, r2 = _circle_3points_xy([(x1+x0)/2, hc], [1*(x1+x0)/4, hm[1]], [x0, he[1]])
        h3, k3, r3 = _circle_3points_xy([(y1+y0)/2, hc], [3*(y1+y0)/4, hm[2]], [y1, he[2]])
        h4, k4, r4 = _circle_3points_xy([(y1+y0)/2, hc], [1*(y1+y0)/4, hm[3]], [y0, he[3]])

    middle = zeros((len(x), 1))

    for i in range(len(x)):
        xi, yi = x[i], y[i]

        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
            # Equation (xi - hx) ** 2 + (hi - kx) ** 2 = rx **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi, ly)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
            # Equation (xi - hy) ** 2 + (hi - ky) ** 2 = ry **2 to find the height of the pointed part (middle of quadrant) with that height one find the equivalent radius
            if he:
                hi = k3 + math.sqrt(r3 ** 2 - (yi - h3) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi, lx)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
            if he:
                hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi, ly)
            if yi <= (y1 + y0)/2:
                zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
            if he:
                hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
            else:
                hi = hc
            ri = _find_r_given_h_l(hi, lx)
            if xi <= (x0 + x1)/2:
                zi = _sqrt((ri)**2 - (xi-(x0+ri))**2)
            else:
                zi = _sqrt((ri)**2 - (xi-(x1-ri))**2)

        else:
            print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))

        middle[i] = zi

    return middle


def pointed_vault_ub_lb_update(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], hc=8.0, he=None, hm=None, tol=1e-6):
    """Update upper and lower bounds of a pointed vault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the arch
    t : float
        Parameter for lower bound in nodes in the boundary
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    he : [float, float, float, float], optional
        Height of the opening mid-span for each of the quadrants, by default None
    hm : [float, float, float, float], optional
        Height of each quadrant center (spadrel), by default None
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    # y1_ub = y1 + thk/2
    # y0_ub = y0 - thk/2
    # x1_ub = x1 + thk/2
    # x0_ub = x0 - thk/2

    y1_lb = y1 - thk/2
    y0_lb = y0 + thk/2
    x1_lb = x1 - thk/2
    x0_lb = x0 + thk/2

    # lx_ub = x1_ub - x0_ub
    # ly_ub = y1_ub - y0_ub

    # lx_lb = x1_lb - x0_lb
    # ly_lb = y1_lb - y0_lb

    lx = x1 - x0
    ly = y1 - y0

    # hc_ub = hc + thk/2
    # hc_lb = hc - thk/2

    # r_x = _find_r_given_h_l(hc, x1 - x0)
    # r_y = _find_r_given_h_l(hc, y1 - y0)

    if he:
        he_ub = he.copy()
        he_lb = he.copy()
        for i in range(len(he)):
            he_ub[i] += thk/2
            he_lb[i] -= thk/2
    if hm:
        raise NotImplementedError()
        hm_ub = hm.copy()
        hm_lb = hm.copy()
        for i in range(len(hm)):
            hm_ub[i] += thk/2
            hm_lb[i] -= thk/2

    if he and hm is None:
        # h1_ub, k1_ub, r1_ub = _circle_3points_xy([x0, he_ub[1]], [(x1+x0)/2, hc_ub], [x1, he_ub[0]])
        # h2_ub, k2_ub, r2_ub = h1_ub, k1_ub, r1_ub
        # h3_ub, k3_ub, r3_ub = _circle_3points_xy([y0, he_ub[3]], [(y1+y0)/2, hc_ub], [y1, he_ub[2]])
        # h4_ub, k4_ub, r4_ub = h3_ub, k3_ub, r3_ub

        # h1_lb, k1_lb, r1_lb = _circle_3points_xy([x0, he_lb[1]], [(x1+x0)/2, hc_lb], [x1, he_lb[0]])
        # h2_lb, k2_lb, r2_lb = h1_lb, k1_lb, r1_lb
        # h3_lb, k3_lb, r3_lb = _circle_3points_xy([y0, he_lb[3]], [(y1+y0)/2, hc_lb], [y1, he_lb[2]])
        # h4_lb, k4_lb, r4_lb = h3_lb, k3_lb, r3_lb

        h1, k1, r1 = _circle_3points_xy([x0, he[1]], [(x1+x0)/2, hc], [x1, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = _circle_3points_xy([y0, he[3]], [(y1+y0)/2, hc], [y1, he[2]])
        h4, k4, r4 = h3, k3, r3

    # elif hm and he:
    #     h1_ub, k1_ub, r1_ub = _circle_3points_xy([(x1+x0)/2, hc_ub], [3*(x1+x0)/4, hm_ub[0]], [x1, he_ub[0]])
    #     h2_ub, k2_ub, r2_ub = _circle_3points_xy([(x1+x0)/2, hc_ub], [1*(x1+x0)/4, hm_ub[1]], [x0, he_ub[1]])
    #     h3_ub, k3_ub, r3_ub = _circle_3points_xy([(y1+y0)/2, hc_ub], [3*(y1+y0)/4, hm_ub[2]], [y1, he_ub[2]])
    #     h4_ub, k4_ub, r4_ub = _circle_3points_xy([(y1+y0)/2, hc_ub], [1*(y1+y0)/4, hm_ub[3]], [y0, he_ub[3]])

    #     h1_lb, k1_lb, r1_lb = _circle_3points_xy([(x1+x0)/2, hc_lb], [3*(x1+x0)/4, hm_lb[0]], [x1, he_lb[0]])
    #     h2_lb, k2_lb, r2_lb = _circle_3points_xy([(x1+x0)/2, hc_lb], [1*(x1+x0)/4, hm_lb[1]], [x0, he_lb[1]])
    #     h3_lb, k3_lb, r3_lb = _circle_3points_xy([(y1+y0)/2, hc_lb], [3*(y1+y0)/4, hm_lb[2]], [y1, he_lb[2]])
    #     h4_lb, k4_lb, r4_lb = _circle_3points_xy([(y1+y0)/2, hc_lb], [1*(y1+y0)/4, hm_lb[3]], [y0, he_lb[3]])

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t

    for i in range(len(x)):
        xi, yi = x[i], y[i]

        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
            if he:
                # hi_ub = k1_ub + math.sqrt(r1_ub ** 2 - (xi - h1_ub) ** 2)  # is in fact hi
                # hi_lb = k1_lb + math.sqrt(r1_lb ** 2 - (xi - h1_lb) ** 2)
                hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, ly_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, ly_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, ly)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_y + thk/2  # This only works for pointed (he only)
            # ri_lb = r_y - thk/2  # This only works for pointed (he only)
            if yi <= (y1 + y0)/2:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
            if he:
                # hi_ub = k3_ub + math.sqrt(r3_ub ** 2 - (yi - h3_ub) ** 2)
                # hi_lb = k3_lb + math.sqrt(r3_lb ** 2 - (yi - h3_lb) ** 2)
                hi = k3 + math.sqrt(r3 ** 2 - (yi - h3) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, lx_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, lx_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, lx)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_x + thk/2  # This only works for pointed (he only)
            # ri_lb = r_x - thk/2  # This only works for pointed (he only)
            if xi <= (x0 + x1)/2:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x1-ri))**2)

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
            if he:
                # hi_ub = k2_ub + math.sqrt(r2_ub ** 2 - (xi - h2_ub) ** 2)
                # hi_lb = k2_lb + math.sqrt(r2_lb ** 2 - (xi - h2_lb) ** 2)
                hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_lb = _find_r_given_h_l(hi_lb, ly_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = _find_r_given_h_l(hi_ub, ly_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, ly)
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_y + thk/2  # This only works for pointed (he only)
            # ri_lb = r_y - thk/2  # This only works for pointed (he only)
            if yi <= (y1 + y0)/2:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y1-ri))**2)

        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
            if he:
                # hi_ub = k4_ub + math.sqrt(r4_ub ** 2 - (yi - h4_ub) ** 2)
                # hi_lb = k4_lb + math.sqrt(r4_lb ** 2 - (yi - h4_lb) ** 2)
                hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, lx_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, lx_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, lx)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_x + thk/2  # This only works for pointed (he only)
            # ri_lb = r_x - thk/2  # This only works for pointed (he only)
            if xi <= (x0 + x1)/2:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x1-ri))**2)
        else:
            print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))

        if ((yi) > y1_lb and ((xi) > x1_lb or (xi) < x0_lb)) or ((yi) < y0_lb and ((xi) > x1_lb or (xi) < x0_lb)):
            lb[i] = - 1*t

    return ub, lb


def pointed_vault_dub_dlb(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], hc=8.0, he=None, hm=None, tol=1e-6):
    """Computes the sensitivities of upper and lower bounds in the x, y coordinates and thickness specified.

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the arch
    t : float
        Parameter for lower bound in nodes in the boundary
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    he : [float, float, float, float], optional
        Height of the opening mid-span for each of the quadrants, by default None
    hm : [float, float, float, float], optional
        Height of each quadrant center (spadrel), by default None
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    dub : array
        Values of the sensitivities for the upper bound in the points
    dlb : array
        Values of the sensitivities for the lower bound in the points
    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    # y1_ub = y1 + thk/2
    # y0_ub = y0 - thk/2
    # x1_ub = x1 + thk/2
    # x0_ub = x0 - thk/2

    y1_lb = y1 - thk/2
    y0_lb = y0 + thk/2
    x1_lb = x1 - thk/2
    x0_lb = x0 + thk/2

    # lx_ub = x1_ub - x0_ub
    # ly_ub = y1_ub - y0_ub

    # lx_lb = x1_lb - x0_lb
    # ly_lb = y1_lb - y0_lb

    lx = x1 - x0
    ly = y1 - y0

    # hc_ub = hc + thk/2
    # hc_lb = hc - thk/2

    # r_x = _find_r_given_h_l(hc, x1 - x0)
    # r_y = _find_r_given_h_l(hc, y1 - y0)

    if he:
        he_ub = he.copy()
        he_lb = he.copy()
        for i in range(len(he)):
            he_ub[i] += thk/2
            he_lb[i] -= thk/2
    if hm:
        raise NotImplementedError()
        hm_ub = hm.copy()
        hm_lb = hm.copy()
        for i in range(len(hm)):
            hm_ub[i] += thk/2
            hm_lb[i] -= thk/2

    if he and hm is None:
        # h1_ub, k1_ub, r1_ub = _circle_3points_xy([x0, he_ub[1]], [(x1+x0)/2, hc_ub], [x1, he_ub[0]])
        # h2_ub, k2_ub, r2_ub = h1_ub, k1_ub, r1_ub
        # h3_ub, k3_ub, r3_ub = _circle_3points_xy([y0, he_ub[3]], [(y1+y0)/2, hc_ub], [y1, he_ub[2]])
        # h4_ub, k4_ub, r4_ub = h3_ub, k3_ub, r3_ub

        # h1_lb, k1_lb, r1_lb = _circle_3points_xy([x0, he_lb[1]], [(x1+x0)/2, hc_lb], [x1, he_lb[0]])
        # h2_lb, k2_lb, r2_lb = h1_lb, k1_lb, r1_lb
        # h3_lb, k3_lb, r3_lb = _circle_3points_xy([y0, he_lb[3]], [(y1+y0)/2, hc_lb], [y1, he_lb[2]])
        # h4_lb, k4_lb, r4_lb = h3_lb, k3_lb, r3_lb

        h1, k1, r1 = _circle_3points_xy([x0, he[1]], [(x1+x0)/2, hc], [x1, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = _circle_3points_xy([y0, he[3]], [(y1+y0)/2, hc], [y1, he[2]])
        h4, k4, r4 = h3, k3, r3

    # elif hm and he:
    #     h1_ub, k1_ub, r1_ub = _circle_3points_xy([(x1+x0)/2, hc_ub], [3*(x1+x0)/4, hm_ub[0]], [x1, he_ub[0]])
    #     h2_ub, k2_ub, r2_ub = _circle_3points_xy([(x1+x0)/2, hc_ub], [1*(x1+x0)/4, hm_ub[1]], [x0, he_ub[1]])
    #     h3_ub, k3_ub, r3_ub = _circle_3points_xy([(y1+y0)/2, hc_ub], [3*(y1+y0)/4, hm_ub[2]], [y1, he_ub[2]])
    #     h4_ub, k4_ub, r4_ub = _circle_3points_xy([(y1+y0)/2, hc_ub], [1*(y1+y0)/4, hm_ub[3]], [y0, he_ub[3]])

    #     h1_lb, k1_lb, r1_lb = _circle_3points_xy([(x1+x0)/2, hc_lb], [3*(x1+x0)/4, hm_lb[0]], [x1, he_lb[0]])
    #     h2_lb, k2_lb, r2_lb = _circle_3points_xy([(x1+x0)/2, hc_lb], [1*(x1+x0)/4, hm_lb[1]], [x0, he_lb[1]])
    #     h3_lb, k3_lb, r3_lb = _circle_3points_xy([(y1+y0)/2, hc_lb], [3*(y1+y0)/4, hm_lb[2]], [y1, he_lb[2]])
    #     h4_lb, k4_lb, r4_lb = _circle_3points_xy([(y1+y0)/2, hc_lb], [1*(y1+y0)/4, hm_lb[3]], [y0, he_lb[3]])

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t
    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))

    for i in range(len(x)):
        xi, yi = x[i], y[i]

        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
            if he:
                # hi_ub = k1_ub + math.sqrt(r1_ub ** 2 - (xi - h1_ub) ** 2)  # is in fact hi
                # hi_lb = k1_lb + math.sqrt(r1_lb ** 2 - (xi - h1_lb) ** 2)
                hi = k1 + math.sqrt(r1 ** 2 - (xi - h1) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, ly_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, ly_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, ly)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_y + thk/2  # This only works for pointed (he only)
            # ri_lb = r_y - thk/2  # This only works for pointed (he only)
            if yi <= (y1 + y0)/2:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y1-ri))**2)
            dub[i] = 1/2 * ri_ub/ub[i]
            dlb[i] = - 1/2 * ri_lb/lb[i]

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
            if he:
                # hi_ub = k3_ub + math.sqrt(r3_ub ** 2 - (yi - h3_ub) ** 2)
                # hi_lb = k3_lb + math.sqrt(r3_lb ** 2 - (yi - h3_lb) ** 2)
                hi = k3 + math.sqrt(r3 ** 2 - (yi - h3) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, lx_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, lx_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, lx)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_x + thk/2  # This only works for pointed (he only)
            # ri_lb = r_x - thk/2  # This only works for pointed (he only)
            if xi <= (x0 + x1)/2:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x1-ri))**2)
            dub[i] = 1/2 * ri_ub/ub[i]
            dlb[i] = - 1/2 * ri_lb/lb[i]

        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
            if he:
                # hi_ub = k2_ub + math.sqrt(r2_ub ** 2 - (xi - h2_ub) ** 2)
                # hi_lb = k2_lb + math.sqrt(r2_lb ** 2 - (xi - h2_lb) ** 2)
                hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_lb = _find_r_given_h_l(hi_lb, ly_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = _find_r_given_h_l(hi_ub, ly_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, ly)
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_y + thk/2  # This only works for pointed (he only)
            # ri_lb = r_y - thk/2  # This only works for pointed (he only)
            if yi <= (y1 + y0)/2:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (yi-(y1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (yi-(y1-ri))**2)
            dub[i] = 1/2 * ri_ub/ub[i]
            dlb[i] = - 1/2 * ri_lb/lb[i]

        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
            if he:
                # hi_ub = k4_ub + math.sqrt(r4_ub ** 2 - (yi - h4_ub) ** 2)
                # hi_lb = k4_lb + math.sqrt(r4_lb ** 2 - (yi - h4_lb) ** 2)
                hi = k4 + math.sqrt(r4 ** 2 - (yi - h4) ** 2)
            else:
                hi = hc
                # hi_ub = hc_ub
                # hi_lb = hc_lb
            # ri_ub = _find_r_given_h_l(hi_ub, lx_ub)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_lb = _find_r_given_h_l(hi_lb, lx_lb)  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri = _find_r_given_h_l(hi, lx)
            ri_ub = ri + thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            ri_lb = ri - thk/2  # This in the equation ri ** 2 =  (xi - xc_) ** 2 + (zi - zc_) ** 2  -> zc = 0.0 and xc_ = (x0 + x1)/2
            # ri_ub = r_x + thk/2  # This only works for pointed (he only)
            # ri_lb = r_x - thk/2  # This only works for pointed (he only)
            if xi <= (x0 + x1)/2:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x0+ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x0+ri))**2)
            else:
                ub[i] = _sqrt((ri_ub)**2 - (xi-(x1-ri))**2)
                lb[i] = _sqrt((ri_lb)**2 - (xi-(x1-ri))**2)
            dub[i] = 1/2 * ri_ub/ub[i]
            dlb[i] = - 1/2 * ri_lb/lb[i]
        else:
            print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))

        if ((yi) > y1_lb and ((xi) > x1_lb or (xi) < x0_lb)) or ((yi) < y0_lb and ((xi) > x1_lb or (xi) < x0_lb)):
            lb[i] = - 1*t
            dlb[i] = 0.0

    return dub, dlb  # ub, lb


def _find_r_given_h_l(h, length):

    r = h**2/length + length/4

    return r


def _circle_3points_xy(p1, p2, p3):

    x1 = p1[0]
    z1 = p1[1]
    x2 = p2[0]
    z2 = p2[1]
    x3 = p3[0]
    z3 = p3[1]

    x12 = x1 - x2
    x13 = x1 - x3
    z12 = z1 - z2
    z13 = z1 - z3
    z31 = z3 - z1
    z21 = z2 - z1
    x31 = x3 - x1
    x21 = x2 - x1

    sx13 = x1**2 - x3**2
    sz13 = z1**2 - z3**2
    sx21 = x2**2 - x1**2
    sz21 = z2**2 - z1**2

    f = ((sx13) * (x12) + (sz13) * (x12) + (sx21) * (x13) + (sz21) * (x13)) / (2 * ((z31) * (x12) - (z21) * (x13)))
    g = ((sx13) * (z12) + (sz13) * (z12) + (sx21) * (z13) + (sz21) * (z13)) / (2 * ((x31) * (z12) - (x21) * (z13)))
    c = - x1 ** 2 - z1 ** 2 - 2 * g * x1 - 2 * f * z1
    h = - g
    k = - f
    r2 = h * h + k * k - c
    r = math.sqrt(r2)

    # print('h: ', h, 'k: ', k, 'r: ', r)

    return h, k, r


def _sqrt(x):
    try:
        sqrt_x = math.sqrt(x)
    except BaseException:
        if x > -10e4:
            sqrt_x = math.sqrt(abs(x))
        else:
            sqrt_x = 0.0
            print('Problems to sqrt: ', x)
    return sqrt_x
