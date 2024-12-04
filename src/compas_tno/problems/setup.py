from math import isnan
from typing import TYPE_CHECKING

from numpy import append
from numpy import array
from numpy import diag
from numpy import vstack
from numpy import zeros

from compas_tno.algorithms import apply_sag
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.algorithms import q_from_variables
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.problems import adapt_problem_to_sym_and_fixed_diagram
from compas_tno.problems import adapt_problem_to_sym_diagram
from compas_tno.problems import callback_create_json
from compas_tno.problems import callback_save_json
from compas_tno.problems import constr_wrapper
from compas_tno.problems import initialise_problem_general
from compas_tno.problems import initialize_loadpath
from compas_tno.problems import initialize_tna
from compas_tno.problems import objective_selector
from compas_tno.problems import sensitivities_wrapper
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import compute_edge_stiffness
from compas_tno.utilities import compute_form_initial_lengths
from compas_tno.utilities import set_b_constraint

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis


def set_up_general_optimisation(analysis: "Analysis"):
    """Set up a nonlinear optimisation problem.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object set up for optimise.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    shape = analysis.shape

    printout = optimiser.settings.get("printout", True)
    plot = optimiser.settings.get("plot", False)
    axis_symmetry = optimiser.settings.get("axis_sym", None)
    sym_loads = optimiser.settings.get("sym_loads", False)
    fjac = optimiser.settings.get("jacobian", False)
    starting_point = optimiser.settings.get("starting_point", "current")
    find_inds = optimiser.settings.get("find_inds", False)
    tol_inds = optimiser.settings.get("tol_inds", None)
    method_ind = optimiser.settings.get("method_ind", "QR")
    qmin = optimiser.settings.get("qmin", -1e4)
    qmax = optimiser.settings.get("qmax", +1e-8)
    features = optimiser.settings.get("features", [])
    save_iterations = optimiser.settings.get("save_iterations", False)
    solver_convex = optimiser.settings.get("solver_convex", "MATLAB")
    autodiff = optimiser.settings.get("autodiff", False)

    pattern_center = form.parameters.get("center", None)

    if shape:
        thk = shape.datashape.get("thk", None)
        form.attributes["thk"] = thk  # for safety storing the thk. If optimisation is minthk, this will be overwritten
    else:
        thk = 0.50

    objective = optimiser.settings["objective"]
    variables = optimiser.settings["variables"]
    constraints = optimiser.settings["constraints"]

    i_k = form.index_vertex()

    qmin_applied = form.edges_attribute("qmin")
    for qmin_applied_i in qmin_applied:
        if qmin_applied_i is None:
            print("Appied qmin / qmax:", qmin, qmax)
            apply_bounds_on_q(form, qmin=qmin, qmax=qmax)
            break

    problem = optimiser.problem
    if not problem:
        problem = initialise_problem_general(form)

    problem.variables = variables
    problem.constraints = constraints
    problem.features = features
    problem.shape = shape
    problem.thk = thk

    if "update-loads" in features:
        F, V0, V1, V2 = form.tributary_matrices(sparse=False)
    else:
        F, V0, V1, V2 = 4 * [None]

    problem.F = F
    problem.V0 = V0
    problem.V1 = V1
    problem.V2 = V2
    problem.ro = shape.ro

    if starting_point == "current":
        pass

    elif starting_point == "sag":
        apply_sag(form, boundary_force=50.0)  # the issue here is that after the sag the problem.x0, problem.y0 are not updated
        initialize_tna(form)

    elif starting_point == "loadpath":
        printout_loadpath = False  # this need to be a proper verbose level
        initialize_loadpath(form, problem=problem, find_inds=find_inds, solver_convex=solver_convex, printout=printout_loadpath)

    elif starting_point == "relax":
        equilibrium_fdm(form)
        initialize_tna(form)

    elif starting_point == "tna" or starting_point == "TNA":
        initialize_tna(form)

    else:
        print("Warning: define starting point")

    problem.q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)

    if "fixed" in features and "sym" in features:
        # print('\n-------- Initialisation with fixed and sym form --------')
        adapt_problem_to_sym_and_fixed_diagram(
            problem,
            form,
            method=method_ind,
            list_axis_symmetry=axis_symmetry,
            center=pattern_center,
            correct_loads=sym_loads,
            printout=printout,
            tol=tol_inds,
        )

    elif "sym" in features:
        # print('\n-------- Initialisation with sym form --------')
        adapt_problem_to_sym_diagram(
            problem,
            form,
            list_axis_symmetry=axis_symmetry,
            center=pattern_center,
            correct_loads=sym_loads,
            printout=printout,
        )

    elif "fixed" in features:
        # print('\n-------- Initialisation with fixed form --------')
        adapt_problem_to_fixed_diagram(problem, form, method=method_ind, printout=printout, tol=tol_inds)

    else:
        # print('\n-------- Initialisation with no-fixed and no-sym form --------')
        pass

    # If objective is the complementary energy

    if "Ecomp" in objective.split("-"):
        dXb = array(optimiser.settings["support_displacement"])
        problem.dXb = dXb

    if objective == "Ecomp-nonlinear":
        Ecomp_method = optimiser.settings.get("Ecomp_method", "simplified")
        problem.Ecomp_method = Ecomp_method

        if Ecomp_method == "simplified":
            stiff = zeros((problem.m, 1))
            lengths = compute_form_initial_lengths(form)
            k = compute_edge_stiffness(form, lengths=lengths)
            for index, edge in enumerate(form.edges()):
                stiff[index] = 1 / 2 * 1 / k[index] * lengths[index] ** 2
            problem.stiff = stiff

        elif Ecomp_method == "complete":
            raise NotImplementedError()

    # Set specific constraints

    if "reac_bounds" in constraints:
        problem.b = set_b_constraint(form, printout)
    else:
        problem.b = None

    # Select Objetive, Gradient, Costraints and Jacobian.

    fobj, fgrad = objective_selector(objective)

    fconstr = constr_wrapper
    if fjac:
        fjac = sensitivities_wrapper

    # Alternative for autodiff

    if autodiff:
        raise NotImplementedError("Autodifferentiation is currently not available")

    # Select starting point (x0) and max/min for variables

    x0 = problem.q[problem.ind]
    bounds = [[qmin.item(), qmax.item()] for i, (qmin, qmax) in enumerate(zip(problem.qmin, problem.qmax)) if i in problem.ind]

    qid = problem.q[problem.ind]
    problem.q = q_from_variables(qid, problem.B, problem.d)
    problem.X[problem.free] = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)

    error = sum((problem.E.dot(problem.q) - problem.ph) ** 2)
    if error > 0.001:
        print("Warning: Error equilibrium:", error)

    # bounds = [[qmin, qmax]] * problem.k

    if "xyb" in variables:
        xyb0 = problem.X[problem.fixed, :2].flatten("F").reshape((-1, 1))
        x0 = append(x0, xyb0).reshape(-1, 1)
        bounds_x = []
        bounds_y = []
        for i in problem.fixed:
            bounds_x.append([form.vertex_attribute(i_k[i], "xmin"), form.vertex_attribute(i_k[i], "xmax")])
            bounds_y.append([form.vertex_attribute(i_k[i], "ymin"), form.vertex_attribute(i_k[i], "ymax")])
        bounds = bounds + bounds_x + bounds_y

    if "zb" in variables:
        zb0 = problem.X[problem.fixed, 2].flatten("F").reshape((-1, 1))
        x0 = append(x0, zb0).reshape(-1, 1)
        bounds_z = []
        for i in problem.fixed:
            bounds_z.append([form.vertex_attribute(i_k[i], "lb"), form.vertex_attribute(i_k[i], "ub")])
        bounds = bounds + bounds_z

    if "t" in variables:
        min_thk = optimiser.settings.get("min_thk", 0.001)
        max_thk = optimiser.settings.get("max_thk", thk)
        x0 = append(x0, thk).reshape(-1, 1)
        bounds = bounds + [[min_thk, max_thk]]

    if "n" in variables:
        thk0_approx = thk  # shape.datashape['thk']
        print("Thickness approximate:", thk0_approx)
        x0 = append(x0, 0.0).reshape(-1, 1)
        min_limit = 0.0  # /2  # 0.0
        bounds = bounds + [[min_limit, thk0_approx / 2]]
        # bounds = bounds + [[min_limit, thk0_approx/2]]

    if "tub" in variables:
        problem.tub = zeros((problem.n, 1))
        tubmax = form.vertices_attribute("tubmax")
        problem.tubmax = array(tubmax).reshape(problem.n, 1)
        problem.tubmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tub)
        bounds = bounds + list(zip(problem.tubmin, problem.tubmax))

    if "tlb" in variables:
        problem.tlb = zeros((problem.n, 1))
        tlbmax = form.vertices_attribute("tlbmax")
        problem.tlbmax = array(tlbmax).reshape(problem.n, 1)
        problem.tlbmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tlb)
        bounds = bounds + list(zip(problem.tlbmin, problem.tlbmax))

    if "tub_reac" in variables:
        tub_reac = []
        for key in form.vertices_where({"is_fixed": True}):
            tub_reac.append(form.vertex_attribute(key, "tub_reacmax"))
        tub_reac = array(tub_reac)
        tub_reac = vstack([tub_reac[:, 0].reshape(-1, 1), tub_reac[:, 1].reshape(-1, 1)])
        problem.tub_reac = abs(tub_reac)
        problem.tub_reac_min = zeros((2 * problem.nb, 1))
        x0 = append(x0, problem.tub_reac_min)
        bounds = bounds + list(zip(problem.tub_reac_min, problem.tub_reac))

    if "lambdh" in variables:
        lambd0 = 1.0
        direction = optimiser.settings.get("load_direction", None).reshape(-1, 1)
        problem.px0 = direction[: problem.n].reshape(-1, 1)
        problem.py0 = direction[problem.n : 2 * problem.n].reshape(-1, 1)
        max_lambd = optimiser.settings.get("max_lambd", lambd0 * 10)
        min_lambd = 0.0
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    if "lambdv" in variables:
        direction = array(optimiser.settings.get("load_direction", None)).reshape(-1, 1)
        max_lambd = optimiser.settings.get("max_lambd", 100.0)
        min_lambd = 0.0
        lambd0 = 1.0
        problem.pzv = direction
        problem.pz0 = array(form.vertices_attribute("pz")).reshape(-1, 1)
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    if "delta" in variables:
        dX = optimiser.settings["displ_map"]
        max_delta = optimiser.settings.get("max_delta", 1.0)
        min_delta = 0.0
        delta0 = optimiser.settings.get("delta0", 0.0)
        problem.dX = dX
        x0 = append(x0, delta0).reshape(-1, 1)
        bounds = bounds + [[min_delta, max_delta]]

        Ud = diag(problem.C @ dX[:, 0])
        Vd = diag(problem.C @ dX[:, 1])
        Edx = problem.Cit @ Ud
        Edy = problem.Cit @ Vd
        problem.Ed = vstack([Edx, Edy])

    if save_iterations:
        callback_create_json()
        optimiser.callback = callback_save_json
        callback_save_json(x0)  # save staring point to file

    # if plot:
    #     plotter = TNOPlotter(form)
    #     plotter.draw_form_independents()
    #     plotter.show()
    #     if "sym" in features:
    #         plotter = TNOPlotter(form)
    #         plotter.draw_form_sym(print_sym=True)
    #         plotter.show()

    if any([isnan(problem.ub[i]) for i in range(len(problem.ub))]) or any([isnan(problem.lb[i]) for i in range(len(problem.lb))]):
        print("Is Nan for the bounds. Optimisation can not proceed")
        raise ValueError("Check bounds that constraint nodes")

    f0 = fobj(x0, problem)
    g0 = fconstr(x0, problem)

    if fgrad:
        grad = fgrad(x0, problem)
    if fjac:
        jac = fjac(x0, problem)

    for i, key in enumerate(form.vertices()):
        form.vertex_attribute(key, "x", problem.X[i, 0])
        form.vertex_attribute(key, "y", problem.X[i, 1])
        form.vertex_attribute(key, "z", problem.X[i, 2])

    if plot:
        from compas_tno.viewers import Viewer

        view = Viewer(form)
        view.draw_thrust()
        view.show()

    if printout:
        print("-" * 20)
        print("NPL (Non Linear Problem) Data:")
        print("Number of force variables:", len(problem.ind))
        print("Number of variables:", len(x0))
        print("Number of constraints:", len(g0))
        if "funicular" in problem.constraints:
            print("# constraints funicular:", 2 * len(problem.q))
        if "envelopexy" in problem.constraints:
            print("# constraints envelope xy:", 4 * len(problem.X))
        if "envelope" in problem.constraints:
            print("# constraints envelope z:", 2 * len(problem.X))
        if "reac_bounds" in problem.constraints:
            print("# constraints reac_bounds:", 2 * len(problem.fixed))
        if fgrad:
            print("Shape of gradient:", grad.shape)
        if fjac:
            print("Shape of jacobian:", jac.shape)
        print("Init. Objective Value: {0}".format(f0))
        print("Init. Constraints Extremes: {0:.3f} to {1:.3f}".format(max(g0), min(g0)))
        violated = []
        for i in range(len(g0)):
            if g0[i] < 0:
                violated.append(i)
        if violated:
            print("Constraints Violated #:", violated)

    optimiser.fobj = fobj
    optimiser.fconstr = fconstr
    optimiser.fgrad = fgrad
    optimiser.fjac = fjac
    optimiser.x0 = x0
    optimiser.bounds = bounds
    optimiser.f0 = f0
    optimiser.g0 = g0

    optimiser.problem = problem

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis


def set_up_convex_optimisation(analysis: "Analysis"):
    """Set up a convex optimisation problem.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object set up for optimise.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    qmax = optimiser.settings["qmax"]
    qmin = optimiser.settings["qmin"]

    objective = optimiser.settings["objective"]
    variables = optimiser.settings["variables"]
    constraints = optimiser.settings["constraints"]

    # Select Objective

    if objective not in ["loadpath", "feasibility"]:
        print("Warning: Non-convex problem for the objective: ", objective, ". Try changing the objective to 'loadpath' or 'fesibility'.")

    if variables == ["q"]:
        pass
    else:
        print("Warning: Non-convex problem for the variables: ", variables, ". Considering only 'q' instead and assuming coplanar supports (zb=0).")

    if constraints == ["funicular"]:
        pass
    else:
        print("Warning: Non-convex problem for the constraints: ", constraints, ". Considering only 'funicular' instead.")

    apply_bounds_on_q(form, qmin=qmin, qmax=qmax)

    problem = initialise_problem_general(form)
    problem.variables = variables
    problem.constraints = constraints

    optimiser.problem = problem

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis
