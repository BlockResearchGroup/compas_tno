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
from compas_tno.problems import FixedProblem
from compas_tno.problems import FixedSymmetricProblem
from compas_tno.problems import Problem
from compas_tno.problems import SymmetricProblem
from compas_tno.problems import callback_create_json
from compas_tno.problems import callback_save_json
from compas_tno.problems import constr_wrapper
from compas_tno.problems import objective_selector
from compas_tno.problems import sensitivities_wrapper
from compas_tno.utilities import compute_edge_stiffness
from compas_tno.utilities import compute_form_initial_lengths

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis


# =============================================================================
# Helper Functions
# =============================================================================


def _apply_bounds_to_edges(form, qmin, qmax):
    """Apply force density bounds to all edges."""
    if qmin is not None:
        for edge in form.edges_where({"_is_edge": True}):
            form.edge_attribute(edge, "qmin", qmin)
    if qmax is not None:
        for edge in form.edges_where({"_is_edge": True}):
            form.edge_attribute(edge, "qmax", qmax)


def _adapt_problem_to_features(form, problem, features, settings, printout):
    """Adapt the problem to the features."""
    method_ind = settings.get("method_ind", "QR")
    axis_symmetry = settings.get("axis_sym", None)
    sym_loads = settings.get("sym_loads", False)
    tol_inds = settings.get("tol_inds", None)
    pattern_center = form.centroid()

    if "fixed" in features and "sym" in features:
        return FixedSymmetricProblem.from_formdiagram(
            form,
            base_problem=problem,
            method=method_ind,
            list_axis_symmetry=axis_symmetry,
            center=pattern_center,
            correct_loads=sym_loads,
            printout=printout,
            tol=tol_inds,
        )
    elif "sym" in features:
        return SymmetricProblem.from_formdiagram(
            form,
            base_problem=problem,
            list_axis_symmetry=axis_symmetry,
            center=pattern_center,
            correct_loads=sym_loads,
            printout=printout,
        )
    elif "fixed" in features:
        return FixedProblem.from_formdiagram(
            form,
            base_problem=problem,
            method=method_ind,
            printout=printout,
            tol=tol_inds,
        )
    else:
        return problem


def _setup_problem_metadata(problem, envelope, form, variables, constraints, features, thk):
    """Setup basic problem metadata."""
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
    problem.rho = envelope.rho


def _apply_starting_point(form, problem, starting_point, settings):
    """Apply starting point strategy to form diagram."""
    # Local imports to avoid circular dependency
    from compas_tno.solvers.startingpoint import startingpoint_loadpath
    from compas_tno.solvers.startingpoint import startingpoint_sag
    from compas_tno.solvers.startingpoint import startingpoint_tna

    if starting_point == "current":
        pass
    elif starting_point == "sag":
        boundary_force = 50.0
        startingpoint_sag(form, boundary_force=boundary_force)
    elif starting_point == "loadpath":
        printout_loadpath = False
        find_inds = settings.get("find_inds", False)
        solver_convex = settings.get("solver_convex", "CLARABEL")
        startingpoint_loadpath(form, problem=problem, find_inds=find_inds, solver_convex=solver_convex, printout=printout_loadpath)
    elif starting_point == "relax":
        equilibrium_fdm(form)
        startingpoint_tna(form)
    elif starting_point == "tna" or starting_point == "TNA":
        startingpoint_tna(form)
    else:
        print("Warning: define starting point")

    # Update problem.q from form
    problem.q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)


def _setup_objective_specific_params(problem, objective, settings, form):
    """Setup objective-specific parameters."""
    if "Ecomp" in objective.split("-"):
        dXb = array(settings["support_displacement"])
        problem.dXb = dXb

    if objective == "Ecomp-nonlinear":
        Ecomp_method = settings.get("Ecomp_method", "simplified")
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


def _setup_constraint_params(problem, constraints, form):
    """Setup constraint-specific parameters."""
    if "reac_bounds" in constraints:
        fixed = list(form.vertices_where({"is_support": True}))
        array(form.vertices_attribute("b", keys=fixed))
        problem.b = array(form.vertices_attribute("b", keys=fixed))
    else:
        problem.b = None


def _build_variable_vector(problem, form, variables, settings, i_k, thk):
    """Build optimization variable vector (x0) and bounds from variable list."""
    x0 = problem.q[problem.ind]
    bounds = [[qmin_i.item(), qmax_i.item()] for i, (qmin_i, qmax_i) in enumerate(zip(problem.qmin, problem.qmax)) if i in problem.ind]

    # Support position variables (xyb)
    if "xyb" in variables:
        xyb0 = problem.X[problem.fixed, :2].flatten("F").reshape((-1, 1))
        x0 = append(x0, xyb0).reshape(-1, 1)
        bounds_x = []
        bounds_y = []
        for i in problem.fixed:
            bounds_x.append([form.vertex_attribute(i_k[i], "xmin"), form.vertex_attribute(i_k[i], "xmax")])
            bounds_y.append([form.vertex_attribute(i_k[i], "ymin"), form.vertex_attribute(i_k[i], "ymax")])
        bounds = bounds + bounds_x + bounds_y

    # Support height variables (zb)
    if "zb" in variables:
        zb0 = problem.X[problem.fixed, 2].flatten("F").reshape((-1, 1))
        x0 = append(x0, zb0).reshape(-1, 1)
        bounds_z = []
        for i in problem.fixed:
            bounds_z.append([form.vertex_attribute(i_k[i], "lb"), form.vertex_attribute(i_k[i], "ub")])
        bounds = bounds + bounds_z

    # Thickness variable (t)
    if "t" in variables:
        min_thk = settings.get("min_thk", 0.001)
        max_thk = settings.get("max_thk", thk)
        x0 = append(x0, thk).reshape(-1, 1)
        bounds = bounds + [[min_thk, max_thk]]

    # Normal offset variable (n)
    if "n" in variables:
        thk0_approx = thk
        print("Thickness approximate:", thk0_approx)
        x0 = append(x0, 0.0).reshape(-1, 1)
        min_limit = 0.0
        bounds = bounds + [[min_limit, thk0_approx / 2]]

    # Extrados thickness variable (tub)
    if "tub" in variables:
        problem.tub = zeros((problem.n, 1))
        tubmax = form.vertices_attribute("tubmax")
        problem.tubmax = array(tubmax).reshape(problem.n, 1)
        problem.tubmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tub)
        bounds = bounds + list(zip(problem.tubmin, problem.tubmax))

    # Intrados thickness variable (tlb)
    if "tlb" in variables:
        problem.tlb = zeros((problem.n, 1))
        tlbmax = form.vertices_attribute("tlbmax")
        problem.tlbmax = array(tlbmax).reshape(problem.n, 1)
        problem.tlbmin = zeros((problem.n, 1))
        x0 = append(x0, problem.tlb)
        bounds = bounds + list(zip(problem.tlbmin, problem.tlbmax))

    # Reaction thickness variable (tub_reac)
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

    # Horizontal load multiplier (lambdh)
    if "lambdh" in variables:
        lambd0 = 1.0
        direction = settings.get("load_direction", None).reshape(-1, 1)
        problem.px0 = direction[: problem.n].reshape(-1, 1)
        problem.py0 = direction[problem.n : 2 * problem.n].reshape(-1, 1)
        max_lambd = settings.get("max_lambd", lambd0 * 10)
        min_lambd = 0.0
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    # Vertical load multiplier (lambdv)
    if "lambdv" in variables:
        direction = array(settings.get("load_direction", None)).reshape(-1, 1)
        max_lambd = settings.get("max_lambd", 100.0)
        min_lambd = 0.0
        lambd0 = 1.0
        problem.pzv = direction
        problem.pz0 = array(form.vertices_attribute("pz")).reshape(-1, 1)
        x0 = append(x0, lambd0).reshape(-1, 1)
        bounds = bounds + [[min_lambd, max_lambd]]

    # Displacement variable (delta)
    if "delta" in variables:
        dX = settings["displ_map"]
        max_delta = settings.get("max_delta", 1.0)
        min_delta = 0.0
        delta0 = settings.get("delta0", 0.0)
        problem.dX = dX
        x0 = append(x0, delta0).reshape(-1, 1)
        bounds = bounds + [[min_delta, max_delta]]

        Ud = diag(problem.C @ dX[:, 0])
        Vd = diag(problem.C @ dX[:, 1])
        Edx = problem.Cit @ Ud
        Edy = problem.Cit @ Vd
        problem.Ed = vstack([Edx, Edy])

    return x0, bounds


def _compute_initial_values(problem, fobj, fconstr, fgrad, fjac, x0):
    """Compute initial objective and constraint values and validate bounds."""
    f0 = fobj(x0, problem)
    g0 = fconstr(x0, problem)

    grad = None
    jac = None
    if fgrad:
        grad = fgrad(x0, problem)
    if fjac:
        jac = fjac(x0, problem)

    if any([isnan(problem.ub[i]) for i in range(len(problem.ub))]) or any([isnan(problem.lb[i]) for i in range(len(problem.lb))]):
        print("Is Nan for the bounds. Optimisation can not proceed")
        raise ValueError("Check bounds that constraint nodes")

    return f0, g0, grad, jac


# def _update_form_coordinates(form, problem):
#     """Update form diagram with current problem coordinates."""
#     for i, key in enumerate(form.vertices()):
#         form.vertex_attribute(key, "x", problem.X[i, 0])
#         form.vertex_attribute(key, "y", problem.X[i, 1])
#         form.vertex_attribute(key, "z", problem.X[i, 2])


def _print_diagnostics(problem, x0, g0, fgrad, fjac, grad, jac, f0):
    """Print problem setup diagnostics."""
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

    if fgrad and grad is not None:
        print("Shape of Gradient:", grad.shape)
    if fjac and jac is not None:
        print("Shape of Jacobian:", jac.shape)

    print("Initial Objective Value: {0}".format(f0))
    print("Initial Constraints Extremes: {0:.3f} to {1:.3f}".format(max(g0), min(g0)))

    violated = []
    for i in range(len(g0)):
        if g0[i] < 0:
            violated.append(i)
    if violated:
        print("# Constraints Violated at Start:", len(violated))


# =============================================================================
# Main Setup Function
# =============================================================================


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

    # Extract components
    form = analysis.formdiagram
    envelope = analysis.envelope
    optimiser = analysis.optimiser
    settings = optimiser.settings
    problem = optimiser.problem

    # Extract key settings
    printout = settings.get("printout", True)
    features = settings.get("features", [])
    save_iterations = settings.get("save_iterations", False)
    autodiff = settings.get("autodiff", False)
    qmin = settings.get("qmin", -1e4)
    qmax = settings.get("qmax", +1e-8)

    objective = settings["objective"]
    variables = settings["variables"]
    constraints = settings["constraints"]

    thk = envelope.thickness
    i_k = form.index_vertex()

    # 1. Apply force density bounds to edges
    _apply_bounds_to_edges(form, qmin, qmax)

    # 2. Create or get problem instance
    if not problem:
        problem = Problem.from_formdiagram(form)
    
    problem = _adapt_problem_to_features(form, problem, features, settings, printout)

    # 3. Setup problem metadata
    _setup_problem_metadata(problem, envelope, form, variables, constraints, features, thk)

    # 4. Starting point should already be applied before calling this function
    # Update problem.q from current form state
    problem.q = array([form.edge_attribute((u, v), "q") for u, v in form.edges_where({"_is_edge": True})]).reshape(-1, 1)

    # 5. Setup objective-specific parameters
    _setup_objective_specific_params(problem, objective, settings, form)

    # 6. Select objective and gradient functions
    fobj, fgrad = objective_selector(objective)

    # 7. Setup constraint-specific parameters
    _setup_constraint_params(problem, constraints, form)

    # 8. Setup constraint and jacobian functions
    fconstr = constr_wrapper
    fjac_enabled = settings.get("jacobian", False)
    fjac = sensitivities_wrapper if fjac_enabled else None

    # 9. Initialize q variables and update problem geometry (before adding other variables)
    qid = problem.q[problem.ind]
    problem.q = q_from_variables(qid, problem.B, problem.d)
    problem.X[problem.free] = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)

    # 10. Verify equilibrium
    error = sum((problem.E.dot(problem.q) - problem.ph) ** 2)
    if error > 0.001:
        print("Warning: Error equilibrium:", error)

    if autodiff:
        raise NotImplementedError("Autodifferentiation is currently not available")

    # 11. Build full variable vector (x0) and bounds (adds q and other variables)
    x0, bounds = _build_variable_vector(problem, form, variables, settings, i_k, thk)

    # 12. Setup callbacks if needed
    if save_iterations:
        callback_create_json()
        optimiser.callback = callback_save_json
        callback_save_json(x0)

    # 13. Compute initial values
    f0, g0, grad, jac = _compute_initial_values(problem, fobj, fconstr, fgrad, fjac, x0)

    # 14. Print diagnostics
    if printout:
        _print_diagnostics(problem, x0, g0, fgrad, fjac, grad, jac, f0)

    # 15. Store functions and variables in problem
    problem.fobj = fobj
    problem.fconstr = fconstr
    problem.fgrad = fgrad
    problem.fjac = fjac
    problem.x0 = x0
    problem.bounds = bounds
    problem.f0 = f0
    problem.g0 = g0

    # 16. Store problem reference in optimiser
    optimiser.problem = problem

    return analysis


def set_up_base_problem(analysis: "Analysis"):
    """Set up the base problem structure (lightweight, no optimization).

    This creates the fundamental problem structure from the form diagram,
    including equilibrium matrices, connectivity, and bounds. No starting
    point computation or optimization is performed.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : :class:`~compas_tno.analysis.Analysis`
        Analysis object with base problem created in optimiser.problem.

    """

    form = analysis.formdiagram
    optimiser = analysis.optimiser
    qmax = optimiser.settings.get("qmax", 1e-8)
    qmin = optimiser.settings.get("qmin", -1e4)

    variables = optimiser.settings.get("variables", ["q"])
    constraints = optimiser.settings.get("constraints", ["funicular"])

    # 1. Apply force density bounds to edges
    _apply_bounds_to_edges(form, qmin, qmax)

    # 2. Create base Problem instance (no features for base problem)
    problem = Problem.from_formdiagram(form)
    problem.variables = variables
    problem.constraints = constraints

    # 3. Store in optimiser
    optimiser.problem = problem

    return analysis


# Legacy alias for backward compatibility
def set_up_convex_optimisation(analysis: "Analysis"):
    """Legacy alias for set_up_base_problem. Use set_up_base_problem instead."""
    return set_up_base_problem(analysis)
