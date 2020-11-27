from compas_tno.problems.problems import initialise_problem

from compas_tno.problems.objectives import f_min_loadpath
from compas_tno.problems.objectives import f_min_thrust
from compas_tno.problems.objectives import f_max_thrust
from compas_tno.problems.objectives import f_target
from compas_tno.problems.objectives import f_constant
from compas_tno.problems.objectives import f_reduce_thk
from compas_tno.problems.objectives import f_tight_crosssection

from compas_tno.problems.derivatives import gradient_fmin
from compas_tno.problems.derivatives import gradient_fmax
from compas_tno.problems.derivatives import gradient_feasibility
from compas_tno.problems.derivatives import gradient_reduce_thk
from compas_tno.problems.derivatives import gradient_tight_crosssection

from compas_tno.problems.derivatives import sensitivities_wrapper
from compas_tno.problems.derivatives import sensitivities_wrapper_inequalities

from compas_tno.plotters import plot_sym_inds

from compas_tno.problems.constraints import constr_wrapper
from compas_tno.problems.constraints import constr_wrapper_inequalities

from compas.datastructures import mesh_bounding_box_xy

from numpy import append
from numpy import array

__all__ = [
    'set_up_nonlinear_optimisation',
    'set_up_convex_optimisation',
]


def set_up_nonlinear_optimisation(analysis):
    """ Set up a nonlinear optimisation problem.

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object set up for optimise.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    shape = analysis.shape

    indset = form.attributes['indset']
    find_inds = optimiser.data.get('find_inds', True)
    printout = optimiser.data.get('printout', True)
    qmax = optimiser.data.get('qmax', 1e+6)
    qmin = optimiser.data.get('qmin', 0.0)  # This qmin sometimes is 1e-6...
    plot = optimiser.data.get('plot', False)

    objective = optimiser.data['objective']
    variables = optimiser.data['variables']
    constraints = optimiser.data['constraints']

    i_k = form.index_key()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    # Set specific constraints

    if 'reac_bounds' in constraints:
        b = set_b_constraint(form, printout)
    else:
        b = None

    if 'cracks' in constraints:
        cracks_lb, cracks_ub = set_cracks_constraint(form, printout)
    else:
        cracks_lb, cracks_ub = None, None

    if 'rollers' in constraints:
        max_rol_rx, max_rol_ry = set_rollers_constraint(form, printout)
    else:
        max_rol_rx, max_rol_ry = None, None

    if 'joints' in constraints:
        joints = set_joints_constraint(form, printout)
    else:
        joints = None

    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in constraints):
        Asym = set_symmetry_constraint(form, constraints, printout, plot=plot)
    else:
        Asym = None

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind,
            ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin,
            constraints,  max_rol_rx, max_rol_ry, Asym, variables, shape)

    # Select Objetive and Gradient

    if objective == 'loadpath':
        fobj = f_min_loadpath
        fgrad = None
    if objective == 'target' or objective == 'bestfit':
        fobj = f_target
        fgrad = None
    if objective == 'min':
        fobj = f_min_thrust
        fgrad = gradient_fmin
    if objective == 'max':
        fobj = f_max_thrust
        fgrad = gradient_fmax
    if objective == 'feasibility':
        fobj = f_constant
        fgrad = gradient_feasibility
    if objective == 't':  # analytical reduce thickness
        fobj = f_reduce_thk
        fgrad = gradient_reduce_thk
    if objective == 's':  # tight UB and LB 0 -> 1/2
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    if objective == 'n':  # vector n offset the surfaces -> larger the better (higher GSF)
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection

    # Select Constraints and Jacobian (w/ equalities in IPOPT and only Inequalities in SLSQP)

    if optimiser.data['solver'] == 'slsqp' or optimiser.data['solver'] == 'SLSQP':
        fconstr = constr_wrapper_inequalities
        fjac = sensitivities_wrapper_inequalities
    else:
        fconstr = constr_wrapper
        fjac = sensitivities_wrapper

    # Definition of the Variables and starting point

    if 'ind' in variables:
        x0 = q[ind]
        bounds = [[qmin, qmax]] * k
    else:
        x0 = q
        bounds = [[qmin, qmax]] * len(q)

    if 'zb' in variables:
        x0 = append(x0, z[fixed]).reshape(-1, 1)
        zb_bounds = [[form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')] for i in fixed]
        bounds = bounds + zb_bounds

    if 't' in variables:
        min_thk = optimiser.data.get('min_thk', 0.001)
        thk = optimiser.data.get('thk', shape.data['thk'])
        x0 = append(x0, thk).reshape(-1, 1)
        bounds = bounds + [[min_thk, thk]]

    if 's' in variables:
        x0 = append(x0, 0.0).reshape(-1, 1)
        bounds = bounds + [[-1.0, 0.5]]

    if 'n' in variables:
        thk0_approx = min(ub - lb)
        print('Thickness approximate:', thk0_approx)
        x0 = append(x0, 0.0).reshape(-1, 1)
        min_limit = - thk0_approx#/2  # 0.0
        bounds = bounds + [[min_limit, thk0_approx/2]]

    f0 = fobj(x0, *args)
    g0 = fconstr(x0, *args)

    if fgrad:
        grad = fgrad(x0, *args)
    if fjac:
        jac = fjac(x0, *args)

    if printout:
        print('-'*20)
        print('NPL (Non Linear Problem) Data:')
        print('Total of Independents:', len(ind))
        print('Number of Variables:', len(x0))
        print('Number of Constraints:', len(g0))
        if fgrad:
            print('Shape of Gradient:', grad.shape)
        if fjac:
            print('Shape of Jacobian:', jac.shape)
        print('Init. Objective Value: {0}'.format(f0))
        print('Init. Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))
        print('Constraint Mapping: dep: 0-{0} | z: {0}-{1} | reac-bound: {1}-{2}'.format(len(dep), len(dep)+2*len(z), len(dep)+2*len(z)+len(fixed)))
        violated = []
        for i in range(len(g0)):
            if g0[i] < 0:
                violated.append(i)
        if violated:
            print('Constraints Violated #:', violated)

    optimiser.fobj = fobj
    optimiser.fconstr = fconstr
    optimiser.fgrad = fgrad
    optimiser.fjac = fjac
    optimiser.args = args
    optimiser.x0 = x0
    optimiser.bounds = bounds
    optimiser.f0 = f0
    optimiser.g0 = g0

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis


def set_b_constraint(form, printout):
    b = []
    for key in form.vertices_where({'is_fixed': True}):
        try:
            [b_] = form.vertex_attributes(key, 'b')
            b.append(b_)
        except:
            pass
    b = array(b)
    if printout:
        print('Reaction bounds active in : {0} joints'.format(len(b)))
    return b


def set_joints_constraint(form, printout):
    try:
        joints = form.attributes['joints']
    except:
        joints = None
    if printout and joints:
        print('Constraints on the Joints set for {0} contacts.'.format(len(joints)))
    return joints


def set_cracks_constraint(form, printout):
    try:
        cracks_lb, cracks_ub = form.attributes['cracks']
        if printout:
            print('Cracks Definition')
            print(cracks_lb, cracks_ub)
    except:
        cracks_lb = []
        cracks_ub = []
    if printout:
        print('Constraints on cracks activated in {0} lb and {1} ub.'.format(len(cracks_lb), len(cracks_ub)))
    return cracks_lb, cracks_ub


def set_rollers_constraint(form, printout):
    max_rol_rx = []
    max_rol_ry = []
    for key in form.vertices_where({'rol_x': True}):
        max_rol_rx.append(form.vertex_attribute(key, 'max_rx'))
    for key in form.vertices_where({'rol_y': True}):
        max_rol_ry.append(form.vertex_attribute(key, 'max_ry'))
    if printout:
        print('Constraints on rollers activated in {0} in x and {1} in y.'.format(len(max_rol_rx), len(max_rol_ry)))
    return array(max_rol_rx).reshape(-1, 1), array(max_rol_ry).reshape(-1, 1)


def set_symmetry_constraint(form, constraints, printout, plot=False):

    horizontal_only = False
    vertical_only = False
    corners = mesh_bounding_box_xy(form)
    xs = [point[0] for point in corners]
    ys = [point[1] for point in corners]
    xc = (max(xs) - min(xs))/2
    yc = (max(ys) - min(ys))/2
    if 'symmetry-horizontal' in constraints:
        horizontal_only = True
    if 'symmetry-vertical' in constraints:
        vertical_only = True
    form.apply_symmetry(center=[xc, yc, 0.0], horizontal_only=horizontal_only, vertical_only=vertical_only)
    Asym = form.assemble_symmetry_matrix()
    if printout:
        print('Calculated and found symmetry from point:', xc, yc)
        print('Resulted in Asym Matrix Shape:', Asym.shape)
        print('Unique independents:', form.number_of_sym_independents())
    if plot:
        plot_sym_inds(form).show()

    return Asym


def set_up_convex_optimisation(analysis):
    """ Set up a nonlinear optimisation problem.

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object set up for optimise.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    indset = form.attributes['indset']
    find_inds = optimiser.data['find_inds']
    printout = optimiser.data['printout']
    qmax = optimiser.data['qmax']

    k_i = form.key_index()
    i_uv = form.index_uv()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y = args[:31]

    # Problem specifities

    args_cvx = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i)
    objective = optimiser.data['objective']

    # Select Objective

    if objective == 'loadpath' or objective == 'feasibility':
        pass
    else:
        print('Error! Non-covex problem for the objective: ', objective, '. Try changing the objective to \'loadpath\' or \'fesibility\'.')

    # Definition of the Variables and starting point

    variables = optimiser.data['variables']

    if 'ind' in variables and 'zb' not in variables:
        pass
    else:
        print('Error! Non-covex problem for the variables: ', variables, '. Try allow for only \'ind\' variables.')

    constraints = optimiser.data['constraints']

    if constraints == ['funicular']:
        pass
    else:
        print('Error! Non-covex problem for the constraints: ', constraints, '. Try allow for only \'funicular\' variables.')

    optimiser.args = args_cvx

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis
