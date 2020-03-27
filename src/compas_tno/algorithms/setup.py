from compas_tno.algorithms.problems import initialise_problem

from compas_tno.algorithms.objectives import f_min_loadpath
from compas_tno.algorithms.objectives import f_min_thrust
from compas_tno.algorithms.objectives import f_max_thrust
from compas_tno.algorithms.objectives import f_target
from compas_tno.algorithms.objectives import f_constant

from compas_tno.algorithms.constraints import f_compression
from compas_tno.algorithms.constraints import f_ub_lb
from compas_tno.algorithms.constraints import f_joints
from compas_tno.algorithms.constraints import f_cracks
from compas_tno.algorithms.constraints import constr_wrapper

from numpy import append
from numpy import array
from numpy import asscalar

from compas_tno.diagrams import FormDiagram

__all__ =[
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
    indset = form.attributes['indset']
    find_inds = optimiser.data['find_inds']
    printout = optimiser.data['printout']
    objective = optimiser.data['objective']
    variables = optimiser.data['variables']
    qmax = optimiser.data['qmax']
    qmin = optimiser.data['qmin']
    constraints = optimiser.data['constraints']

    i_k = form.index_key()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    # Set constraints

    if 'reac_bounds' in constraints:
        b = set_b_constraint(form, True, True)
    else:
        b = None

    if 'cracks' in constraints:
        cracks_lb, cracks_ub = set_cracks_constraint(form, True, True)
    else:
        cracks_lb, cracks_ub = None, None

    if 'joints' in constraints:
        joints = set_joints_constraint(form, True)
    else:
        joints = None

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind,
            ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints)

    fconstr = constr_wrapper
    # fconstr = f_ub_lb

    # Select Objetive

    if objective == 'loadpath':
        fobj = f_min_loadpath
    if objective == 'target':
        fobj = f_target
    if objective == 'min':
        fobj = f_min_thrust
    if objective == 'max':
        fobj = f_max_thrust
    if objective == 'feasibility':
        fobj = f_constant

    # Definition of the Variables and starting point

    if 'ind' in variables and 'zb' in variables:
        x0 = q[ind]
        zb_bounds = [[form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')] for i in fixed]
        bounds = [[-10e-6, qmax]] * k + zb_bounds
        x0 = append(x0, z[fixed]).reshape(-1, 1)
    else:
        x0 = q[ind]
        bounds = [[-10e-6, qmax]] * k

    print('Total of Independents:', len(ind))
    print('Number of Variables:', len(x0))
    f0 = fobj(x0, *args)
    g0 = fconstr(x0, *args)
    print('Number of Constraints:', len(g0))

    print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
    print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))

    optimiser.fobj = fobj
    optimiser.fconstr = fconstr
    optimiser.args = args
    optimiser.x0 = x0
    optimiser.bounds = bounds
    optimiser.f0 = f0
    optimiser.g0 = g0

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis


def set_b_constraint(form, bmax, printout):
    if bmax:
        b = []
        for key in form.vertices_where({'is_fixed': True}):
            try:
                [b_] = form.vertex_attributes(key, 'b')
                b.append(b_)
            except:
                pass
        b = array(b)
    else:
        b = None
    if printout and bmax:
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


def set_cracks_constraint(form, cracks, printout):
    if cracks:
        try:
            cracks_lb, cracks_ub = form.attributes['cracks']
            if printout:
                print('Cracks Definition')
                print(cracks_lb, cracks_ub)
        except:
            cracks_lb = []
            cracks_ub = []
    else:
        cracks_lb = []
        cracks_ub = []
    print('Constraints on cracks activated in {0} lb and {1} ub.'.format(len(cracks_lb), len(cracks_ub)))
    return cracks_lb, cracks_ub


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


def set_up_nonlinear_optimisation_tensor(analysis):
    # Work in Progress, will unify the calculations with tensor...
    # """ Set up a nonlinear optimisation problem.

    # Parameters
    # ----------
    # obj : analysis
    #     Analysis object with information about optimiser, form and shape.

    # Returns
    # -------
    # obj : analysis
    #     Analysis object set up for optimise.

    # """

    # form = analysis.form
    # optimiser = analysis.optimiser
    # indset = form.attributes['indset']
    # find_inds = optimiser.data['find_inds']
    # printout = optimiser.data['printout']
    # objective = optimiser.data['objective']
    # variables = optimiser.data['variables']
    # qmax = optimiser.data['qmax']
    # qmin = optimiser.data['qmin']
    # constraints = optimiser.data['constraints']

    # i_k = form.index_key()

    # args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
    # q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    # # Set constraints

    # if 'reac_bounds' in constraints:
    #     b = set_b_constraint(form, True, True)
    # else:
    #     b = None

    # if 'cracks' in constraints:
    #     cracks_lb, cracks_ub = set_cracks_constraint(form, True, True)
    # else:
    #     cracks_lb, cracks_ub = None, None

    # if 'joints' in constraints:
    #     joints = set_joints_constraint(form, True)
    # else:
    #     joints = None

    # args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind,
    #         ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints)

    # fconstr = constr_wrapper
    # # fconstr = f_ub_lb

    # # Select Objetive

    # if objective == 'loadpath':
    #     fobj = f_min_loadpath
    # if objective == 'target':
    #     fobj = f_target
    # if objective == 'min':
    #     fobj = f_min_thrust
    # if objective == 'max':
    #     fobj = f_max_thrust
    # if objective == 'feasibility':
    #     fobj = f_constant

    # # Definition of the Variables and starting point

    # if 'ind' in variables and 'zb' in variables:
    #     x0 = q[ind]
    #     zb_bounds = [[form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')] for i in fixed]
    #     bounds = [[-10e-6, qmax]] * k + zb_bounds
    #     x0 = append(x0, z[fixed]).reshape(-1, 1)
    # else:
    #     x0 = q[ind]
    #     bounds = [[-10e-6, qmax]] * k

    # print('Total of Independents:', len(ind))
    # print('Number of Variables:', len(x0))
    # f0 = fobj(x0, *args)
    # g0 = fconstr(x0, *args)
    # print('Number of Constraints:', len(g0))

    # print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
    # print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))

    # optimiser.fobj = fobj
    # optimiser.fconstr = fconstr
    # optimiser.args = args
    # optimiser.x0 = x0
    # optimiser.bounds = bounds
    # optimiser.f0 = f0
    # optimiser.g0 = g0

    # analysis.form = form
    # analysis.optimiser = optimiser

    return analysis
