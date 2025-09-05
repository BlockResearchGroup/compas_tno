from math import isnan
from typing import TYPE_CHECKING

from numpy import append
from numpy import array
from numpy import diag
from numpy import vstack
from numpy import zeros

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
from compas_tno.problems import objective_selector
from compas_tno.problems import sensitivities_wrapper
from compas_tno.problems import startingpoint_loadpath
from compas_tno.problems import startingpoint_sag
from compas_tno.problems import startingpoint_tna
from compas_tno.utilities import compute_edge_stiffness
from compas_tno.utilities import compute_form_initial_lengths

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

    # =============================================================================
    # Initialisation
    # =============================================================================

    form = analysis.formdiagram
    envelope = analysis.envelope
    optimiser = analysis.optimiser

    printout = optimiser.settings.get("printout", True)
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
    solver_convex = optimiser.settings.get("solver_convex", "CLARABEL")
    autodiff = optimiser.settings.get("autodiff", False)

    pattern_center = form.centroid()

    thk = analysis.envelope.thickness

    objective = optimiser.settings["objective"]
    variables = optimiser.settings["variables"]
    constraints = optimiser.settings["constraints"]

    i_k = form.index_vertex()

    if qmin is not None:
        for edge in form.edges_where({"_is_edge": True}):
            form.edge_attribute(edge, "qmin", qmin)
    if qmax is not None:
        for edge in form.edges_where({"_is_edge": True}):
            form.edge_attribute(edge, "qmax", qmax)

    problem = optimiser.problem
    if not problem:
        problem = initialise_problem_general(form)

    problem.variables = variables
    problem.constraints = constraints
    problem.features = features
    problem.envelope = envelope
    problem.thk = thk

    if "update-loads" in features:
        F, V0, V1, V2 = form.tributary_matrices(sparse=False)
    else:
        F, V0, V1, V2 = 4 * [None]

    problem.F = F
    problem.V0 = V0
    problem.V1 = V1
    problem.V2 = V2
    problem.rho = analysis.envelope.rho

    # =============================================================================
    # Starting point
    # =============================================================================

    if starting_point == "current":
        pass

    elif starting_point == "sag":
        boundary_force = 50.0
        startingpoint_sag(form, boundary_force=boundary_force)

    elif starting_point == "loadpath":
        printout_loadpath = False  # this need to be a proper verbose level
        startingpoint_loadpath(form, problem=problem, find_inds=find_inds, solver_convex=solver_convex, printout=printout_loadpath)

    elif starting_point == "relax":
        equilibrium_fdm(form)
        startingpoint_tna(form)

    elif starting_point == "tna" or starting_point == "TNA":
        startingpoint_tna(form)

    else:
        print("Warning: define starting point")

    problem.q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)

    # =============================================================================
    # Adaptation to symmetry and fixed form
    # =============================================================================

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

    # =============================================================================
    # Objective
    # =============================================================================

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

    # Select Objetive, Gradient, Costraints and Jacobian.

    fobj, fgrad = objective_selector(objective)

    # =============================================================================
    # Constraints
    # =============================================================================

    if "reac_bounds" in constraints:
        # double check if reac_bounds have been applied
        fixed = list(form.vertices_where({"is_support": True}))
        array(form.vertices_attribute("b", keys=fixed))
        problem.b = array(form.vertices_attribute("b", keys=fixed))
    else:
        problem.b = None

    fconstr = constr_wrapper
    if fjac:
        fjac = sensitivities_wrapper

    # =============================================================================
    # Starting point
    # =============================================================================

    x0 = problem.q[problem.ind]
    bounds = [[qmin.item(), qmax.item()] for i, (qmin, qmax) in enumerate(zip(problem.qmin, problem.qmax)) if i in problem.ind]

    problem.q = q_from_variables(x0, problem.B, problem.d)
    problem.X[problem.free] = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)

    error = sum((problem.E.dot(problem.q) - problem.ph) ** 2)
    if error > 0.001:
        print("Warning: Error equilibrium:", error)

    if autodiff:
        raise NotImplementedError("Autodifferentiation is currently not available")

    # =============================================================================
    # Variables
    # =============================================================================

    # If the supports can move in xy, add them as variables
    if "xyb" in variables:
        xyb0 = problem.X[problem.fixed, :2].flatten("F").reshape((-1, 1))
        x0 = append(x0, xyb0).reshape(-1, 1)
        bounds_x = []
        bounds_y = []
        for i in problem.fixed:
            bounds_x.append([form.vertex_attribute(i_k[i], "xmin"), form.vertex_attribute(i_k[i], "xmax")])
            bounds_y.append([form.vertex_attribute(i_k[i], "ymin"), form.vertex_attribute(i_k[i], "ymax")])
        bounds = bounds + bounds_x + bounds_y

    # If the supports can move in z, add it as a variable (x0) and add the bounds
    if "zb" in variables:
        zb0 = problem.X[problem.fixed, 2].flatten("F").reshape((-1, 1))
        x0 = append(x0, zb0).reshape(-1, 1)
        bounds_z = []
        for i in problem.fixed:
            bounds_z.append([form.vertex_attribute(i_k[i], "lb"), form.vertex_attribute(i_k[i], "ub")])
        bounds = bounds + bounds_z

    # If the thickness can change, add it as a variable (x0) and add the bounds
    if "t" in variables:
        min_thk = optimiser.settings.get("min_thk", 0.001)
        max_thk = optimiser.settings.get("max_thk", thk)
        x0 = append(x0, thk).reshape(-1, 1)
        bounds = bounds + [[min_thk, max_thk]]

    # If the thickness can change, add it as a variable (x0) and add the bounds
    if "n" in variables:
        thk0_approx = thk  # shape.parameters['thk']
        print("Thickness approximate:", thk0_approx)
        x0 = append(x0, 0.0).reshape(-1, 1)
        min_limit = 0.0  # /2  # 0.0
        bounds = bounds + [[min_limit, thk0_approx / 2]]
        # bounds = bounds + [[min_limit, thk0_approx/2]]

    # If a variable to increase the 'z' of the extrados is considered
    if "tub" in variables:
        problem.tub = zeros((problem.n, 1))
        tubmax = form.vertices_attribute("tubmax")
        problem.tubmax = array(tubmax).reshape(problem.n, 1)
        problem.tubmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tub)
        bounds = bounds + list(zip(problem.tubmin, problem.tubmax))

    # If a variable to decrease the 'z' of the intrados is considered
    if "tlb" in variables:
        problem.tlb = zeros((problem.n, 1))
        tlbmax = form.vertices_attribute("tlbmax")
        problem.tlbmax = array(tlbmax).reshape(problem.n, 1)
        problem.tlbmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tlb)
        bounds = bounds + list(zip(problem.tlbmin, problem.tlbmax))

    # If a variable to increase the 'z' of the reaction forces is considered
    if "tub_reac" in variables:
        tub_reac = []
        for key in form.vertices_where({"is_support": True}):
            tub_reac.append(form.vertex_attribute(key, "tub_reacmax"))
        tub_reac = array(tub_reac)
        tub_reac = vstack([tub_reac[:, 0].reshape(-1, 1), tub_reac[:, 1].reshape(-1, 1)])
        problem.tub_reac = abs(tub_reac)
        problem.tub_reac_min = zeros((2 * problem.nb, 1))
        x0 = append(x0, problem.tub_reac_min)
        bounds = bounds + list(zip(problem.tub_reac_min, problem.tub_reac))

    # If a variable to increase horizontal load multiplier 'lambdh' is considered
    if "lambdh" in variables:
        lambd0 = 1.0
        direction = optimiser.settings.get("load_direction", None).reshape(-1, 1)
        problem.px0 = direction[: problem.n].reshape(-1, 1)
        problem.py0 = direction[problem.n : 2 * problem.n].reshape(-1, 1)
        max_lambd = optimiser.settings.get("max_lambd", lambd0 * 10)
        min_lambd = 0.0
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    # If a variable to increase vertical load multiplier 'lambdv' of the applied loads is considered
    if "lambdv" in variables:
        direction = array(optimiser.settings.get("load_direction", None)).reshape(-1, 1)
        max_lambd = optimiser.settings.get("max_lambd", 100.0)
        min_lambd = 0.0
        lambd0 = 1.0
        problem.pzv = direction
        problem.pz0 = array(form.vertices_attribute("pz")).reshape(-1, 1)
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    # If a variable to displace the nodes horizontally is considered
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

    if printout:
        print("-" * 20)
        print("TNO v.2.0")
        print("Non Linear Problem Data:")
        print("Number of q variables:", len(problem.ind))
        print("Number of total variables:", len(x0))
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
            print("Shape of Gradient:", grad.shape)
        if fjac:
            print("Shape of Jacobian:", jac.shape)
        print("Initial Objective Value: {0}".format(f0))
        print("Initial Constraints Extremes: {0:.3f} to {1:.3f}".format(max(g0), min(g0)))
        violated = []
        for i in range(len(g0)):
            if g0[i] < 0:
                violated.append(i)
        if violated:
            print("# Constraints Violated at Start:", len(violated))

    optimiser.fobj = fobj
    optimiser.fconstr = fconstr
    optimiser.fgrad = fgrad
    optimiser.fjac = fjac
    optimiser.x0 = x0
    optimiser.bounds = bounds
    optimiser.f0 = f0
    optimiser.g0 = g0
    optimiser.problem = problem

    # analysis.form = form  # No need I think, it's updated in place
    # analysis.optimiser = optimiser

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

    form = analysis.formdiagram
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

    # Apply bounds on the edges' force densities (apply_bounds_on_q)
    if isinstance(qmin, list):
        for i, edge in enumerate(form.edges_where({"_is_edge": True})):
            form.edge_attribute(edge, "qmin", qmin[i])
            form.edge_attribute(edge, "qmax", qmax[i])
    else:
        for i, edge in enumerate(form.edges_where({"_is_edge": True})):
            form.edge_attribute(edge, "qmin", qmin)
            form.edge_attribute(edge, "qmax", qmax)

    problem = initialise_problem_general(form)
    problem.variables = variables
    problem.constraints = constraints
    # problem.thk = None

    # analysis.optimiser.fconstr = None

    optimiser.problem = problem
    # analysis.form = form
    # analysis.optimiser = optimiser

    return analysis
