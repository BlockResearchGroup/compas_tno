import math
from numpy import array


__all__ = [
    'check_envelope_constraints',
    'distance_target',
    'rectangular_smoothing_constraints',
    'assign_cracks',
    'rollers_on_openings',
    'set_b_constraint',
    'set_rollers_constraint'
]


def check_envelope_constraints(form, tol=1e-6):
    """Check if the envelope constraints are respected and return a ``penalty > 0`` if not.

    Parameters
    ----------
    form : FormDiagram
        The form diagram with the envelope as nodal attributes.
    tol : float, optional
        Tolerance for verifying the constraints, by default ``1e-6``.

    Returns
    -------
    penalty: float
        The penalty computed as the sum square of the exceeding value at each node.
    """
    penalty = 0.0

    for key in form.vertices():
        _, _, z = form.vertex_coordinates(key)
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        if z < lb - tol:
            outside = (lb - z)
            penalty += outside**2
        if z > ub + tol:
            outside = (z - ub)
            penalty += outside**2

    print('The penalty in the constraints is {0:.3f}'.format(penalty))

    return penalty


def distance_target(form):
    """Returns a measure of the distance among the TN and the target surface.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to check.

    Returns
    -------
    distance: float
        The square distance is computed.

    """

    dist = 0

    for key in form.vertices():
        _, _, z = form.vertex_coordinates(key)
        targ = form.vertex_attribute(key, 'target')
        dist += (z - targ) ** 2

    return dist


def rectangular_smoothing_constraints(form, xy_span=[[0, 10], [0, 10]]):
    """Constraint for smoothing the form diagram in a rectangular boundary.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to check.
    xy_span: list, optional
        The rectangular footprint to smooth, by default [[0, 10], [0, 10]]

    Returns
    -------
    distance: float
        The square distance is computed.

    """

    [[x0, x1], [y0, y1]] = xy_span
    cons = {key: None for key in form.vertices()}
    line_top = [[x0, y1, 0], [x1, y1, 0]]
    line_bottom = [[x0, y0, 0], [x1, y0, 0]]
    line_left = [[x0, y0, 0], [x0, y1, 0]]
    line_right = [[x1, y0, 0], [x1, y1, 0]]
    for key in form.vertices_on_boundary():
        x, y, z = form.vertex_coordinates(key)
        if x == x0:
            cons[key] = line_left
        elif x == x1:
            cons[key] = line_right
        elif y == y0:
            cons[key] = line_bottom
        elif y == y1:
            cons[key] = line_top
    for key in form.vertices_where({'is_fixed': True}):
        cons[key] = form.vertex_coordinates(key)
    return cons


def assign_cracks(form, dx=[[0.50, 0.55]], dy=[[-0.1, 0.1]], type=['top']):
    """Assign cracks on a form diagram to the nodes desired

    Parameters
    ----------
    form : FormDiagram
        ForceDiagram to constraint
    dx : list, optional
        List of the range on the x-coordinates of the nodes to be constrained, by default [[0.50, 0.55]]
    dy : list, optional
        List of the range on the y-coordinates of the nodes to be constrained, by default [[-0.1, 0.1]]
    type : list, optional
        List with the type of constraint to applye to the vertices (top or bottom), by default ['top']

    Returns
    -------
    form : FormDiagram
        ForceDiagram with the constraints in attribute 'cracks'.
    """

    k_i = form.key_index()
    cracks_ub = []
    cracks_lb = []

    # count = 0

    for key in form.vertices():
        x, y, _ = form.vertex_coordinates(key)
        for i in range(len(dx)):
            if dx[i][0] <= x <= dx[i][1] and dy[i][0] <= y <= dy[i][1]:
                if type[i] == 'top':
                    cracks_ub.append(k_i[key])
                if type[i] == 'bottom':
                    cracks_lb.append(k_i[key])

    form.attributes['cracks'] = (cracks_lb, cracks_ub)

    return form


def rollers_on_openings(form, xy_span=[[0.0, 10.0], [0.0, 10.0]], max_f=5.0, constraint_directions='all'):
    """ Apply rollers to the rectangular boundaries of the pattern.

    Parameters
    ----------
    form : FormDiagram
        ForceDiagram to constraint.
    xy_span: list, optional
        The rectangular footprint to smooth, by default [[0, 10], [0, 10]]
    max_f : float, optional
        The maximum force allowed in each roller, by default 5.0
    constraint_directions : str
        Define the restriction in ``all`` directions or only ``x``or only ``y``.

    Returns
    -------
    form : FormDiagram
        ForceDiagram with the rollers assigned.
    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    bndr = form.vertices_on_boundary()

    for key in bndr:
        if form.vertex_attribute(key, 'is_fixed') is False:
            x, y, _ = form.vertex_coordinates(key)
            if x == x1 and (constraint_directions in ['all', 'x']):
                form.vertex_attribute(key, 'rol_x', True)
                form.vertex_attribute(key, 'max_rx', max_f)
            if x == x0 and (constraint_directions in ['all', 'x']):
                form.vertex_attribute(key, 'rol_x', True)
                form.vertex_attribute(key, 'max_rx', max_f)
            if y == y1 and (constraint_directions in ['all', 'y']):
                form.vertex_attribute(key, 'rol_y', True)
                form.vertex_attribute(key, 'max_ry', max_f)
            if y == y0 and (constraint_directions in ['all', 'y']):
                form.vertex_attribute(key, 'rol_y', True)
                form.vertex_attribute(key, 'max_ry', max_f)

    return form


def circular_joints(form, x0=None, xf=None, blocks=18, thk=0.5, t=0.0, tol=1e-3):
    """Deprecated function to assign constraints in the joints among blocks instead of in the z's.

    Parameters
    ----------
    form : FormDiagram
        The form diagram
    x0 : float, optional
        The initial ordinate of the arch, by default None
    xf : float, optional
        The final ordinate of the arch, by default None
    blocks : int, optional
        Number of blocks, by default 18
    thk : float, optional
        The thickness of the arch, by default 0.5
    t : float, optional
        Lower bound in the vertices of the boundary, by default 0.0
    tol : float, optional
        The tolerance, by default 1e-3

    Returns
    -------
    form : FormDiagram
        The form diagram with the joints
    """

    k_i = form.key_index()

    if x0 is None or xf is None:
        x = []
        for key in form.vertices():
            x.append(form.vertex_coordinates(key)[0])
        x0 = min(x)
        xf = max(x)
    y = 0.0

    xc = (xf+x0)/2
    r = xf - xc
    ri = r - thk/2
    re = r + thk/2
    form.attributes['Re'] = re
    form.attributes['Ri'] = ri
    print('SpanMid: {0:.2} m / SpanInt: {1:.2} m / SpanExt: {2:.2} m / Thickness: {3:.4} m / Ratio t/Ri: {4:.4} m / Ratio t/R: {5:.4} m / Number of Blocks: {6}'.format(2 *
          r, 2*ri, 2*re, thk, (thk/ri), (thk/r), blocks))

    njoints = blocks+1
    joints = {}
    for j in range(njoints):
        theta = j/blocks*math.pi
        xi = xc + ri * math.cos(theta)  # takeout
        zi = ri * math.sin(theta) + tol
        xe = xc + re * math.cos(theta)
        ze = re * math.sin(theta) + tol  # take out
        xmax = max(xi, xe)
        xmin = min(xi, xe)
        possible_edges = []
        for (u, v) in form.edges():
            xu, xv = form.vertex_coordinates(u)[0], form.vertex_coordinates(v)[0]
            if max(xu, xv) >= xmin and min(xu, xv) <= xmax:
                possible_edges.append(tuple(sorted([k_i[u], k_i[v]])))
                if form.vertex_attribute(u, 'is_fixed'):
                    possible_edges.append(tuple(sorted([-k_i[u], k_i[u]])))
                if form.vertex_attribute(v, 'is_fixed'):
                    possible_edges.append(tuple(sorted([-k_i[v], k_i[v]])))
        joints[j] = [[xi, y, zi], [xe, y, ze], set(possible_edges)]
        print(joints[j])
    form.attributes['joints'] = joints

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        zt = math.sqrt(r**2 - (x-xc)**2)
        ze = math.sqrt(re**2 - (x-xc)**2) - t
        form.vertex_attribute(key, 'target', value=zt)
        zi2 = ri**2 - (x-xc)**2
        if zi2 < 0:
            zi = 0 - t
        else:
            zi = math.sqrt(zi2) - t
        if form.vertex_attribute(key, 'is_fixed'):
            form.vertex_attribute(key, 'lb', value=None)
            form.vertex_attribute(key, 'ub', value=None)
        else:
            form.vertex_attribute(key, 'lb', value=zi)
            form.vertex_attribute(key, 'ub', value=ze)
        # form.vertex_attribute(key,'z',value=ze)
        if form.vertex_attribute(key, 'is_fixed'):
            form.vertex_attribute(key, 'b', value=[thk/2, 0.0])
        if x == x0:
            form.attributes['tmax'] = ze

    return form


def set_b_constraint(form, printout=False):
    """Helper to set the ``b`` constraint on the reactions of the structure.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to apply the constraints
    printout : bool, optional
        If prints are added to thr screen, by default False

    Returns
    -------
    array (nb x 2)
        The limits for the reaction force.
    """
    b = []
    for key in form.vertices_where({'is_fixed': True}):
        try:
            [b_] = form.vertex_attributes(key, 'b')
            b.append(b_)
        except BaseException:
            pass
    b = array(b)
    if printout:
        print('Reaction bounds active in : {0} joints'.format(len(b)))
    return b


def set_rollers_constraint(form, printout=False):
    """Helper to constraints on rollers.

    Parameters
    ----------
    form : FormDiagram
        The form diagram to apply the constraints
    printout : bool, optional
        If prints are added to thr screen, by default False

    Returns
    -------
    array (nb x 2)
        The limits for the reaction force.
    """
    max_rol_rx = []
    max_rol_ry = []
    for key in form.vertices_where({'rol_x': True}):
        max_rol_rx.append(form.vertex_attribute(key, 'max_rx'))
    for key in form.vertices_where({'rol_y': True}):
        max_rol_ry.append(form.vertex_attribute(key, 'max_ry'))
    if printout:
        print('Constraints on rollers activated in {0} in x and {1} in y.'.format(len(max_rol_rx), len(max_rol_ry)))
    return array(max_rol_rx).reshape(-1, 1), array(max_rol_ry).reshape(-1, 1)
