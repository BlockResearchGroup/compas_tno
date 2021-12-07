from compas_tno.problems import initialise_problem
from compas_tno.problems import initialise_problem_general

from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.problems import adapt_problem_to_sym_diagram
from compas_tno.problems import adapt_problem_to_sym_and_fixed_diagram

from compas_tno.problems import f_min_thrust_general
from compas_tno.problems import f_max_thrust_general
from compas_tno.problems import f_bestfit_general
from compas_tno.problems import f_horprojection_general
from compas_tno.problems import f_loadpath_general
from compas_tno.problems import f_complementary_energy
from compas_tno.problems import f_complementary_energy_nonlinear
from compas_tno.problems import f_max_section

from compas_tno.problems import gradient_fmin_general
from compas_tno.problems import gradient_fmax_general
from compas_tno.problems import gradient_bestfit_general
from compas_tno.problems import gradient_horprojection_general
from compas_tno.problems import gradient_loadpath_general
from compas_tno.problems import gradient_complementary_energy
from compas_tno.problems import gradient_complementary_energy_nonlinear
from compas_tno.problems import gradient_max_section

from compas_tno.problems import f_constant
from compas_tno.problems import f_reduce_thk
from compas_tno.problems import f_tight_crosssection

from compas_tno.problems import callback_save_json
from compas_tno.problems import callback_create_json

from compas_tno.problems import initialize_loadpath
from compas_tno.problems import initialize_tna

from compas_tno.algorithms import apply_sag
from compas_tno.algorithms import z_from_form

from compas_tno.problems import gradient_feasibility
from compas_tno.problems import gradient_reduce_thk
from compas_tno.problems import gradient_tight_crosssection

from compas_tno.problems import sensitivities_wrapper_general

from compas_tno.problems import constr_wrapper_general

from compas_tno.plotters import plot_symmetry
from compas_tno.plotters import plot_symmetry_vertices
from compas_tno.plotters import plot_independents
from compas_tno.plotters import plot_form

from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import compute_form_initial_lengths
from compas_tno.utilities import compute_edge_stiffness
from compas_tno.utilities import compute_average_edge_stiffness

from numpy import append
from numpy import array
from numpy import zeros
from numpy import vstack


__all__ = [
    'set_up_general_optimisation',
    'set_up_convex_optimisation',
]


def set_up_general_optimisation(analysis):
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

    printout = optimiser.settings.get('printout', True)
    plot = optimiser.settings.get('plot', False)
    # thickness_type = optimiser.settings.get('thickness_type', 'constant')
    axis_symmetry = optimiser.settings.get('axis_symmetry', None)
    sym_loads = optimiser.settings.get('sym_loads', False)
    fjac = optimiser.settings.get('jacobian', False)
    starting_point = optimiser.settings.get('starting_point', 'current')
    qmin = optimiser.settings.get('qmin', -1e+4)
    qmax = optimiser.settings.get('qmax', +1e+8)
    features = optimiser.settings.get('features', [])
    save_iterations = optimiser.settings.get('save_iterations', False)
    pattern_center = form.parameters.get('center', None)

    if shape:
        thk = shape.datashape['thk']
    else:
        thk = 0.20

    objective = optimiser.settings['objective']
    variables = optimiser.settings['variables']
    constraints = optimiser.settings['constraints']

    i_k = form.index_key()

    M = initialise_problem_general(form)
    M.variables = variables
    M.constraints = constraints
    M.features = features
    M.shape = shape
    M.thk = thk

    if starting_point == 'current':
        pass
    elif starting_point == 'sag':
        apply_sag(form)
        M.q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)
    elif starting_point == 'loadpath':
        initialize_loadpath(form, problem=M)
        M.q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)
    elif starting_point == 'relax':
        z_from_form(form)
        M.q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)
    elif starting_point == 'tna' or starting_point == 'TNA':
        initialize_tna(form)
    else:
        print('Warning: define starting point')

    if 'fixed' in features and 'sym' in features:
        # print('\n-------- Initialisation with fixed and sym form --------')
        adapt_problem_to_sym_and_fixed_diagram(M, form, list_axis_symmetry=axis_symmetry, center=pattern_center, correct_loads=sym_loads, printout=printout)
    elif 'sym' in features:
        # print('\n-------- Initialisation with sym form --------')
        adapt_problem_to_sym_diagram(M, form, list_axis_symmetry=axis_symmetry, center=pattern_center, correct_loads=sym_loads, printout=printout)
    elif 'fixed' in features:
        # print('\n-------- Initialisation with fixed form --------')
        adapt_problem_to_fixed_diagram(M, form, printout=printout)
    else:
        # print('\n-------- Initialisation with no-fixed and no-sym form --------')
        pass

    if save_iterations:
        callback_create_json()
        optimiser.settings['callback'] = callback_save_json

    apply_bounds_on_q(form, qmin=qmin, qmax=qmax)

    # Specific parameters that depend on the objective:

    if 'Ecomp' in objective.split('-'):
        M.dXb = optimiser.settings['support_displacement']

    if objective == 'Ecomp-nonlinear':
        Ecomp_method = optimiser.settings.get('Ecomp_method', 'simplified')
        M.Ecomp_method = Ecomp_method

        if Ecomp_method == 'simplified':
            stiff = zeros((M.m, 1))
            lengths = compute_form_initial_lengths(form)
            k = compute_edge_stiffness(form, lengths=lengths)
            for index, edge in enumerate(form.edges()):
                stiff[index] = 1 / 2 * 1 / k[index] * lengths[index] ** 2
            M.stiff = stiff
        elif Ecomp_method == 'complete':
            k = compute_average_edge_stiffness(form)
            M.stiff = 1/2 * 1/k

    # Set specific constraints

    if 'reac_bounds' in constraints:
        M.b = set_b_constraint(form, printout)
    else:
        M.b = None

    # if 'cracks' in constraints:
    #     cracks_lb, cracks_ub = set_cracks_constraint(form, printout)
    # else:
    #     cracks_lb, cracks_ub = None, None

    # if 'rollers' in constraints:
    #     max_rol_rx, max_rol_ry = set_rollers_constraint(form, printout)
    # else:
    #     max_rol_rx, max_rol_ry = None, None

    # shape.data['thickness_type'] = thickness_type  # this argument is useful if minimising thickness

    if plot:
        print('Plot of starting point')
        plot_form(form, show_q=False).show()

    # Select Objetive and Gradient

    if objective == 'loadpath':
        fobj = f_loadpath_general
        fgrad = gradient_loadpath_general
    elif objective == 'target' or objective == 'bestfit':
        fobj = f_bestfit_general
        fgrad = gradient_bestfit_general
    elif objective == 'min':
        fobj = f_min_thrust_general
        fgrad = gradient_fmin_general
    elif objective == 'max':
        fobj = f_max_thrust_general
        fgrad = gradient_fmax_general
    elif objective == 'feasibility':
        fobj = f_constant
        fgrad = gradient_feasibility
    elif objective == 'hor_projection':
        fobj = f_horprojection_general
        fgrad = gradient_horprojection_general
    elif objective == 't':  # analytical reduce thickness
        fobj = f_reduce_thk
        fgrad = gradient_reduce_thk
    elif objective == 's':  # tight UB and LB 0 -> 1/2
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    elif objective == 'n':  # vector n offset the surfaces -> larger the better (higher GSF)
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    elif objective == 'lambd':  # vector lambda as hor multiplier larger the better (higher GSF)
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    elif objective == 'Ecomp-linear':  # vector lambda as hor multiplier larger the better (higher GSF)
        fobj = f_complementary_energy
        fgrad = gradient_complementary_energy
    elif objective == 'Ecomp-nonlinear':
        fobj = f_complementary_energy_nonlinear
        fgrad = gradient_complementary_energy_nonlinear
    elif objective == 'max_section':
        fobj = f_max_section
        fgrad = gradient_max_section
    else:
        print('Please, provide a valid objective for the optimisation')
        raise NotImplementedError

    # Select Inequality Constraints and Jacobian

    fconstr = constr_wrapper_general
    if fjac:
        fjac = sensitivities_wrapper_general

    # Select starting point (x0) and max/min for variables

    x0 = M.q[M.ind]
    bounds = [[qmin.item(), qmax.item()] for i, (qmin, qmax) in enumerate(zip(M.qmin, M.qmax)) if i in M.ind]

    # bounds = [[qmin, qmax]] * M.k

    if 'xyb' in variables:
        xyb0 = M.X[M.fixed, :2].flatten('F').reshape((-1, 1))
        x0 = append(x0, xyb0).reshape(-1, 1)
        bounds_x = []
        bounds_y = []
        for i in M.fixed:
            bounds_x.append([form.vertex_attribute(i_k[i], 'xmin'), form.vertex_attribute(i_k[i], 'xmax')])
            bounds_y.append([form.vertex_attribute(i_k[i], 'ymin'), form.vertex_attribute(i_k[i], 'ymax')])
        bounds = bounds + bounds_x + bounds_y

    if 'zb' in variables:
        zb0 = M.X[M.fixed, 2].flatten('F').reshape((-1, 1))
        x0 = append(x0, zb0).reshape(-1, 1)
        bounds_z = []
        for i in M.fixed:
            bounds_z.append([form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')])
        bounds = bounds + bounds_z

    if 't' in variables:
        min_thk = optimiser.settings.get('min_thk', 0.001)
        max_thk = optimiser.settings.get('max_thk', thk)
        print('max thickness', max_thk)
        x0 = append(x0, thk).reshape(-1, 1)
        bounds = bounds + [[min_thk, max_thk]]

    if 'lambd' in variables:
        lambd0 = optimiser.settings.get('lambd', 1.0)
        direction = optimiser.settings.get('lambd-direction', 'x')
        M.px0 = array(form.vertices_attribute('px')).reshape(-1, 1)
        M.py0 = array(form.vertices_attribute('py')).reshape(-1, 1)
        form.apply_horizontal_multiplier(lambd=1.0, direction=direction)
        max_lambd = optimiser.settings.get('max_lambd', lambd0 * 10)
        min_lambd = 0.0
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    # if 's' in variables:
    #     x0 = append(x0, 0.0).reshape(-1, 1)
    #     bounds = bounds + [[-1.0, 0.5]]

    if 'n' in variables:
        thk0_approx = thk  # shape.datashape['thk']
        print('Thickness approximate:', thk0_approx)
        x0 = append(x0, 0.0).reshape(-1, 1)
        min_limit = 0.0  # /2  # 0.0
        bounds = bounds + [[min_limit, thk0_approx/2]]
        # bounds = bounds + [[min_limit, thk0_approx/2]]

    if 'tub' in variables:
        M.tub = zeros((M.n, 1))
        tubmax = form.vertices_attribute('tubmax')
        M.tubmax = array(tubmax).reshape(M.n, 1)
        M.tubmin = zeros((M.n, 1))
        x0 = append(x0, M.tub)
        bounds = bounds + list(zip(M.tubmin, M.tubmax))

    if 'tlb' in variables:
        M.tlb = zeros((M.n, 1))
        tlbmax = form.vertices_attribute('tlbmax')
        M.tlbmax = array(tlbmax).reshape(M.n, 1)
        M.tlbmin = zeros((M.n, 1))
        x0 = append(x0, M.tlb)
        bounds = bounds + list(zip(M.tlbmin, M.tlbmax))

    if 'tub_reac' in variables:
        tub_reac = []
        for key in form.vertices_where({'is_fixed': True}):
            tub_reac.append(form.vertex_attribute(key, 'tub_reacmax'))
        tub_reac = array(tub_reac)
        tub_reac = vstack([tub_reac[:, 0].reshape(-1, 1), tub_reac[:, 1].reshape(-1, 1)])
        M.tub_reac = abs(tub_reac)
        M.tub_reac_min = zeros((2*M.nb, 1))
        x0 = append(x0, M.tub_reac_min)
        bounds = bounds + list(zip(M.tub_reac_min, M.tub_reac))

    # print(M.q)

    # from compas_plotters import MeshPlotter
    # plotter = MeshPlotter(form)
    # plotter.draw_edges(text={edge: round(form.edge_attribute(edge, 'q'), 2) for edge in form.edges()})
    # plotter.show()

    if plot:
        plot_independents(form).show()
        if 'sym' in features:
            plot_symmetry(form).show()
            plot_symmetry_vertices(form).show()

    f0 = fobj(x0, M)
    g0 = fconstr(x0, M)

    # print(max(g0), min(g0))

    if fgrad:
        grad = fgrad(x0, M)
    if fjac:
        jac = fjac(x0, M)

    if printout:
        print('-'*20)
        print('NPL (Non Linear Problem) Data:')
        print('Number of force variables:', len(M.ind))
        print('Number of variables:', len(x0))
        print('Number of constraints:', len(g0))
        if fgrad:
            print('Shape of gradient:', grad.shape)
        if fjac:
            print('Shape of jacobian:', jac.shape)
        print('Init. Objective Value: {0}'.format(f0))
        if objective == 'Ecomp-nonlinear':
            print('Init. Linear Obj Func: {0}'.format(f_complementary_energy(x0, M)))
        print('Init. Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))
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
    optimiser.M = M
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
        except BaseException:
            pass
    b = array(b)
    if printout:
        print('Reaction bounds active in : {0} joints'.format(len(b)))
        print(b)
    return b


def set_joints_constraint(form, printout):
    try:
        joints = form.attributes['joints']
    except BaseException:
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
    except BaseException:
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
    find_inds = optimiser.settings['find_inds']
    printout = optimiser.settings['printout']
    qmax = optimiser.settings['qmax']

    k_i = form.key_index()
    i_uv = form.index_uv()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y = args[:31]

    # Problem specifities

    args_cvx = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i)
    objective = optimiser.settings['objective']

    # Select Objective

    if objective == 'loadpath' or objective == 'feasibility':
        pass
    else:
        print('Error! Non-covex problem for the objective: ', objective, '. Try changing the objective to \'loadpath\' or \'fesibility\'.')

    # Definition of the Variables and starting point

    variables = optimiser.settings['variables']

    if 'ind' in variables and 'zb' not in variables:
        pass
    else:
        print('Error! Non-covex problem for the variables: ', variables, '. Try allow for only \'ind\' variables.')

    constraints = optimiser.settings['constraints']

    if constraints == ['funicular']:
        pass
    else:
        print('Error! Non-covex problem for the constraints: ', constraints, '. Try allow for only \'funicular\' variables.')

    optimiser.args = args_cvx

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis
