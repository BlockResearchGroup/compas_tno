import math
from numpy import arange
from numpy import array
from compas.datastructures import Mesh


# take out ub lb
def pointed_vault_heightfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], hc=8.0, he=None, hm=None, thk=None,  t=0.0, tol=0.00):
    """ Set pointed-vault heights.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].
    hc: float
        Height of the central part of the pointed vault.
    he: list (optional)
        Height of the opening mid-span for each of the quadrants (see Notes).
    hm: list (optional)
        Height of each quadrant mid-span (see Notes).
    ub_lb : bool (optional)
        If True, the thickness will apply and the limits will be stored as attributes 'ub' and 'lb' on the form-diagram
    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface
    tol : float (optional)
        Approximates the equations avoiding negative square-roots.
    set_heights: bool
        If True, the nodes will have the heights 'z' updated to match the pointed arch shape.

    Returns
    -------
    obj
        FormDiagram.

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

    lx = x1 - x0
    ly = y1 - y0
    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0, x1 + lx/density_x, lx/density_x)
    y = arange(y0, y1 + ly/density_y, ly/density_y)

    if he and hm is None:
        h1, k1, r1 = circle_3points_xy([x0, he[1]], [(x1+x0)/2, hc], [x1, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0, he[3]], [(y1+y0)/2, hc], [y1, he[2]])
        h4, k4, r4 = h3, k3, r3
    elif hm and he:
        h1, k1, r1 = circle_3points_xy([(x1+x0)/2, hc], [3*(x1+x0)/4, hm[0]], [x1, he[0]])
        h2, k2, r2 = circle_3points_xy([(x1+x0)/2, hc], [1*(x1+x0)/4, hm[1]], [x0, he[1]])
        h3, k3, r3 = circle_3points_xy([(y1+y0)/2, hc], [3*(y1+y0)/4, hm[2]], [y1, he[2]])
        h4, k4, r4 = circle_3points_xy([(y1+y0)/2, hc], [1*(y1+y0)/4, hm[3]], [y0, he[3]])

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]

            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q1
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

            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q2
                if he:
                    hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi, ly)
                if yi <= (y1 + y0)/2:
                    zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
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

            x1d.append(xi)
            y1d.append(yi)
            z1d.append(zi)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)
            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    middle = Mesh.from_vertices_and_faces(xyz, faces_i)

    # avoid changes in major values of following ub-lb functions
    if he:
        he_ = he.copy()
    else:
        he_ = None
    if hm:
        hm_ = hm.copy()
    else:
        hm_ = None

    extrados = pointed_vault_heightfields_ub(xy_span=xy_span, discretisation=discretisation, hc=hc, he=he_, hm=hm_, thk=thk, tol=tol)

    # avoid changes in major values of following ub-lb functions
    if he:
        he_ = he.copy()
    else:
        he_ = None
    if hm:
        hm_ = hm.copy()
    else:
        hm_ = None

    intrados = pointed_vault_heightfields_lb(xy_span=xy_span, discretisation=discretisation, hc=hc, he=he_, hm=hm_, thk=thk, t=t, tol=tol)

    return intrados, extrados, middle


def pointed_vault_heightfields_ub(xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10,10], hc=8.0, he=None, hm=None, thk=None, tol=0.00):

    y1_init = xy_span[1][1]
    y0_init = xy_span[1][0]
    x1_init = xy_span[0][1]
    x0_init = xy_span[0][0]
    dx = x1_init - x0_init
    dy = y1_init - y0_init
    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0_init, x1_init + dx/density_x, dx/density_x)
    y = arange(y0_init, y1_init + dy/density_y, dy/density_y)

    y1 = xy_span[1][1] + thk / 2
    y0 = xy_span[1][0] - thk / 2
    x1 = xy_span[0][1] + thk / 2
    x0 = xy_span[0][0] - thk / 2

    lx = x1 - x0
    ly = y1 - y0

    hc += thk/2

    if he:
        for i in range(len(he)):
            he[i] += thk/2
    if hm:
        for i in range(len(he)):
            hm[i] += thk/2

    if he and hm is None:
        h1, k1, r1 = circle_3points_xy([x0_init, he[1]], [(x1_init+x0_init)/2, hc], [x1_init, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0_init, he[3]], [(y1_init+y0_init)/2, hc], [y1_init, he[2]])
        h4, k4, r4 = h3, k3, r3
    elif hm and he:
        h1, k1, r1 = circle_3points_xy([(x1_init+x0_init)/2, hc], [3*(x1_init+x0_init)/4, hm[0]], [x1_init, he[0]])
        h2, k2, r2 = circle_3points_xy([(x1_init+x0_init)/2, hc], [1*(x1_init+x0_init)/4, hm[1]], [x0_init, he[1]])
        h3, k3, r3 = circle_3points_xy([(y1_init+y0_init)/2, hc], [3*(y1_init+y0_init)/4, hm[2]], [y1_init, he[2]])
        h4, k4, r4 = circle_3points_xy([(y1_init+y0_init)/2, hc], [1*(y1_init+y0_init)/4, hm[3]], [y0_init, he[3]])

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]

            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q1
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

            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q2
                if he:
                    hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
                else:
                    hi = hc
                ri = _find_r_given_h_l(hi, ly)
                if yi <= (y1 + y0)/2:
                    zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
                else:
                    zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
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

            x1d.append(xi)
            y1d.append(yi)
            z1d.append(zi)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)

            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    extrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return extrados


def pointed_vault_heightfields_lb(xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10,10], hc=8.0, he=None, hm=None, thk=None, t=0.0, tol=0.00):

    y1_init = xy_span[1][1]
    y0_init = xy_span[1][0]
    x1_init = xy_span[0][1]
    x0_init = xy_span[0][0]

    dx = x1_init - x0_init
    dy = y1_init - y0_init

    density_x = discretisation[0]
    density_y = discretisation[1]

    x = arange(x0_init, x1_init + dx/density_x, dx/density_x)
    y = arange(y0_init, y1_init + dy/density_y, dy/density_y)

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    lx = x1 - x0
    ly = y1 - y0

    hc -= thk/2
    if he:
        for i in range(len(he)):
            he[i] -= thk/2
    if hm:
        for i in range(len(hm)):
            hm[i] -= thk/2

    if he and hm:
        h1, k1, r1 = circle_3points_xy([(x1_init+x0_init)/2, hc], [3*(x1_init+x0_init)/4, hm[0]], [x1_init, he[0]])
        h2, k2, r2 = circle_3points_xy([(x1_init+x0_init)/2, hc], [1*(x1_init+x0_init)/4, hm[1]], [x0_init, he[1]])
        h3, k3, r3 = circle_3points_xy([(y1_init+y0_init)/2, hc], [3*(y1_init+y0_init)/4, hm[2]], [y1_init, he[2]])
        h4, k4, r4 = circle_3points_xy([(y1_init+y0_init)/2, hc], [1*(y1_init+y0_init)/4, hm[3]], [y0_init, he[3]])
    elif he and (hm is None):
        h1, k1, r1 = circle_3points_xy([x0_init, he[1]], [(x1_init+x0_init)/2, hc], [x1_init, he[0]])
        h2, k2, r2 = h1, k1, r1
        h3, k3, r3 = circle_3points_xy([y0_init, he[3]], [(y1_init+y0_init)/2, hc], [y1_init, he[2]])
        h4, k4, r4 = h3, k3, r3

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]

            x1d.append(xi)
            y1d.append(yi)

            if ((yi) > y1 and ((xi) > x1 or (xi) < x0)) or ((yi) < y0 and ((xi) > x1 or (xi) < x0)):
                zi = - 1*t
            else:
                if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q1
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

                elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q2
                    if he:
                        hi = k2 + math.sqrt(r2 ** 2 - (xi - h2) ** 2)
                    else:
                        hi = hc
                    ri = _find_r_given_h_l(hi, ly)
                    if yi <= (y1 + y0)/2:
                        zi = _sqrt((ri)**2 - (yi-(y0+ri))**2)
                    else:
                        zi = _sqrt((ri)**2 - (yi-(y1-ri))**2)

                elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
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

            z1d.append(zi)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)

            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    intrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return intrados


def _find_r_given_h_l(h, l):

    r = h**2/l + l/4

    return r


def circle_3points_xy(p1, p2, p3):

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
    except:
        if x > -10e4:
            sqrt_x = math.sqrt(abs(x))
        else:
            sqrt_x = 0.0
            print('Problems to sqrt: ', x)
    return sqrt_x
