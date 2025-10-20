from typing import TYPE_CHECKING
from typing import Optional

import cvxpy as cp
from cvxpy import Minimize
from cvxpy import Problem as CVXProblem
from cvxpy import diag
from cvxpy import matrix_frac

from compas_tno.algorithms import compute_reactions
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import FixedProblem
from compas_tno.problems import Problem as TNOProblem
from compas_tno.solvers.result import SolverResult
from compas_tno.solvers.solver import ConvexSolver

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis


class CVXPYSolver(ConvexSolver):
    """CVXPY solver for convex loadpath optimization problems.

    This solver uses CVXPY to formulate and solve the loadpath minimization
    problem as a convex optimization. It automatically handles both cases:
    - Problems with independent edges (FixedProblem)
    - Problems with all edges as variables (Problem)

    Parameters
    ----------
    settings : dict, optional
        Solver settings. Common options:
        - 'solver_convex' : str - CVXPY backend (default: 'CLARABEL')
          Options: 'CLARABEL', 'MOSEK', 'CVXOPT'
        - 'printout' : bool - Print solver output (default: False)
        - 'update_form' : bool - Update form diagram with solution (default: True)

    Examples
    --------
    >>> solver = CVXPYSolver(settings={"solver_convex": "CLARABEL", "printout": True})
    >>> result = solver.solve(problem)
    >>> if result.success:
    ...     print(f"Loadpath: {result.fopt}")

    Notes
    -----
    This solver is specifically designed for loadpath minimization.
    It formulates the problem as:

    minimize: loadpath = sum(|q_i| * l_i)
    subject to: horizontal equilibrium, force density bounds

    """

    def __init__(self, settings: Optional[dict] = None):
        """Initialize CVXPY solver."""
        super().__init__(settings)

    def _solve_convex(self, problem: "TNOProblem", **kwargs) -> SolverResult:
        """Solve using CVXPY.

        Parameters
        ----------
        problem : Problem
            The problem structure with matrices.
        **kwargs
            Additional options (override settings):
            - solver_convex : str - CVXPY backend
            - printout : bool - Print output
            - update_form : bool - Update FormDiagram

        Returns
        -------
        SolverResult
            The optimization result with optimal force densities.

        """
        solver_convex = kwargs.get("solver_convex", self.settings.get("solver_convex", "CLARABEL"))
        printout = kwargs.get("printout", self.settings.get("printout", False))

        # Determine which formulation to use based on independents
        if len(problem.ind) < problem.m:
            if printout:
                print(f"Calling LP-Optimization via CVXPY with {problem.k} independents")
            fopt, qopt, exitflag, niter, status, sol_time = self._solve_with_independents(problem, solver_convex, printout)
        else:
            if printout:
                print("Calling LP-Optimization via CVXPY with NO independents")
            fopt, qopt, exitflag, niter, status, sol_time = self._solve_without_independents(problem, solver_convex, printout)

        # Check success
        success = status not in ["infeasible", "unbounded"]
        message = f"CVXPY {status}"

        if printout:
            print("\n" + "-" * 50)
            print("LOADPATH OPTIMISATION WITH CVXPY")
            print(f"status    : {status}")
            print(f"fopt (lp) : {fopt:.3f}")
            print(f"n-iter    : {niter}")
            print(f"q range   : {min(qopt):.3f} : {max(qopt):.3f}")
            print(f"sol. time : {sol_time:.3f} sec")
            print("-" * 50 + "\n")

        return SolverResult(
            xopt=qopt.reshape(-1, 1),
            fopt=float(fopt),
            success=success,
            message=message,
            niter=niter,
            time=sol_time,
            exitflag=0 if success else 1,
        )

    def _solve_without_independents(self, problem: "TNOProblem", solver_convex: str, printout: bool):
        """Solve loadpath problem without independent edges.

        Used when all edges are variables (base Problem).

        """
        # Retrieve matrices and vectors
        E = problem.E
        C = problem.C
        Ci = problem.Ci
        Cit = problem.Cit
        Cb = problem.Cb
        pz = problem.P[:, 2].reshape(-1, 1)
        ph = problem.ph
        free = problem.free
        fixed = problem.fixed
        qmin = problem.qmin
        qmax = problem.qmax
        x = problem.x0
        y = problem.y0
        m = problem.m

        # Define CVXPY variables
        q = cp.Variable(m)

        # Objective: minimize loadpath
        fobj = matrix_frac(pz[free], -Cit @ cp.diag(q) @ Ci) - x.T @ C.T @ diag(q) @ Cb @ x[fixed] - y.T @ C.T @ diag(q) @ Cb @ y[fixed]
        objective = Minimize(fobj)

        # Constraints
        horz = E @ q == ph.flatten()
        pos = q >= qmin.flatten()
        maxq = q <= qmax.flatten()
        constraints = [horz, pos, maxq]

        # Solve
        prob = CVXProblem(objective, constraints)
        prob.solve(solver=solver_convex, verbose=printout)

        # Extract results
        fopt = prob.value
        qopt = q.value
        status = prob.status
        niter = prob.solver_stats.num_iters
        sol_time = prob.solver_stats.solve_time
        exitflag = 1 if status not in ["infeasible", "unbounded"] else 0

        return fopt, qopt, exitflag, niter, status, sol_time

    def _solve_with_independents(self, problem: "TNOProblem", solver_convex: str, printout: bool):
        """Solve loadpath problem with independent edges.

        Used for FixedProblem and FixedSymmetricProblem.

        """
        # Retrieve matrices and vectors
        Edinv = problem.Edinv
        Ei = problem.Ei
        ph = problem.ph
        C = problem.C
        Ci = problem.Ci
        Cit = problem.Cit
        Cb = problem.Cb
        pz = problem.P[:, 2].reshape(-1, 1)
        free = problem.free
        fixed = problem.fixed
        dep = problem.dep
        ind = problem.ind
        qmin = problem.qmin
        qmax = problem.qmax
        x = problem.x0
        y = problem.y0
        m = problem.m

        # Define CVXPY variables
        q = cp.Variable(m)

        # Objective: minimize loadpath
        fobj = matrix_frac(pz[free], -Cit @ cp.diag(q) @ Ci) - x.T @ C.T @ diag(q) @ Cb @ x[fixed] - y.T @ C.T @ diag(q) @ Cb @ y[fixed]
        objective = Minimize(fobj)

        # Constraints (including horizontal equilibrium via independents)
        horz = q[dep] == Edinv @ (Ei @ q[ind] - ph.flatten())
        pos = q >= qmin.flatten()
        maxq = q <= qmax.flatten()
        constraints = [horz, pos, maxq]

        # Solve
        prob = CVXProblem(objective, constraints)
        prob.solve(solver=solver_convex, verbose=printout)

        # Extract results
        fopt = prob.value
        qopt = q.value
        status = prob.status
        niter = prob.solver_stats.num_iters
        sol_time = prob.solver_stats.solve_time
        exitflag = 1 if status not in ["infeasible", "unbounded"] else 0

        return fopt, qopt, exitflag, niter, status, sol_time


def run_convex_optimisation(analysis: "Analysis") -> "TNOProblem":
    """Run convex optimisation problem with CVXPY.

    This function is kept for backward compatibility. New code should use CVXPYSolver class.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    result : SolverResult
        The SolverResult object with the key solution information.

    Notes
    -----
    This function wraps the new CVXPYSolver class for backward compatibility.
    For new code, prefer using::

        solver = CVXPYSolver(settings={"solver_convex": "CLARABEL", "printout": True})
        result = solver.solve(problem)

    """
    form = analysis.formdiagram
    optimiser = analysis.optimiser
    problem = optimiser.problem

    # Create solver
    solver = CVXPYSolver(settings=optimiser.settings)

    # Solve
    result = solver.solve(problem)

    # Update problem with optimal solution
    problem.q = result.xopt
    Xfinal = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)
    problem.X[problem.free, 2] = Xfinal[:, 2]

    # Update form diagram
    i_uv = problem.i_uv
    for i, key in enumerate(form.vertices()):
        form.vertex_attribute(key, "z", problem.X[i, 2])

    for c, qi in enumerate(list(problem.q.ravel())):
        edge = i_uv[c]
        li = form.edge_length(edge)
        form.edge_attribute(edge, "q", float(qi))
        form.edge_attribute(edge, "f", float(qi * li))

    compute_reactions(form)

    # # Store results in optimiser
    # optimiser.fopt = result.fopt
    # optimiser.xopt = result.xopt
    # optimiser.exitflag = result.exitflag
    # optimiser.niter = result.niter
    # optimiser.time = result.time
    # optimiser.message = result.message

    return result


def run_loadpath_from_form_CVXPY(form, problem=None, find_inds=False, solver_convex="CLARABEL", printout=False):
    """Run convex optimisation problem with CVXPY directly from the Form Diagram.

    This function is kept for backward compatibility. New code should use CVXPYSolver class.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram object with form containing full information.
    problem : Problem, optional
        The problem with matrices of interest, by default None
    find_inds : bool, optional
        Whether or not independents must be computed before the analysis, by default False
    solver_convex : str, optional
        Solver to use, by default CLARABEL. Options are "CLARABEL", "MOSEK" or "CVXOPT".
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    problem : Problem
        Problem with optimal solution applied.

    Notes
    -----
    For new code, prefer using::

        from compas_tno.problems import FixedProblem
        from compas_tno.solvers import CVXPYSolver

        problem = FixedProblem.from_formdiagram(form, method="QR")
        solver = CVXPYSolver(settings={"solver_convex": "CLARABEL"})
        result = solver.solve(problem)

    """
    # Create appropriate problem type
    if not problem:
        if find_inds:
            problem = FixedProblem.from_formdiagram(form, method="QR")
        else:
            problem = TNOProblem.from_formdiagram(form)

    # Use new CVXPYSolver class
    solver = CVXPYSolver(settings={"solver_convex": solver_convex, "printout": printout})
    result = solver.solve(problem)

    # Update problem and form
    problem.q = result.xopt
    Xfinal = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)
    problem.X[problem.free, 2] = Xfinal[:, 2]

    i_uv = problem.i_uv
    for i, key in enumerate(form.vertices()):
        form.vertex_attribute(key, "z", problem.X[i, 2])

    for c, qi in enumerate(list(problem.q.ravel())):
        edge = i_uv[c]
        li = form.edge_length(edge)
        form.edge_attribute(edge, "q", float(qi))
        form.edge_attribute(edge, "f", float(qi * li))

    compute_reactions(form)

    return result


def call_and_output_CVXPY(form, problem, solver_convex="CLARABEL", printout=False):
    """Call and output the loadpath optimisation with CVXPY.

    This function is kept for backward compatibility. New code should use CVXPYSolver class.

    Parameters
    ----------
    form : FormDiagram
        The form Diagram of the analysis
    problem : Problem
        The Problem with relevant matrices and vectors
    solver_convex : str, optional
        Solver to use, by default CLARABEL. Options are "CLARABEL", "MOSEK" or "CVXOPT".
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    result : SolverResult
        The SolverResult object with the key solution information.

    Notes
    -----
    This is a legacy function. For new code, use CVXPYSolver class directly.

    """
    # Use new CVXPYSolver class
    solver = CVXPYSolver(settings={"solver_convex": solver_convex, "printout": printout})
    result = solver.solve(problem)

    # Update problem and form
    problem.q = result.xopt
    Xfinal = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)
    problem.X[problem.free, 2] = Xfinal[:, 2]

    i_uv = problem.i_uv
    for i, key in enumerate(form.vertices()):
        form.vertex_attribute(key, "z", problem.X[i, 2])

    for c, qi in enumerate(list(problem.q.ravel())):
        edge = i_uv[c]
        li = form.edge_length(edge)
        form.edge_attribute(edge, "q", float(qi))
        form.edge_attribute(edge, "f", float(qi * li))

    compute_reactions(form)

    return result


def call_cvxpy(problem, solver_convex="CLARABEL", printout=False):
    """Legacy helper function - use CVXPYSolver class instead.

    This function is kept for backward compatibility but is no longer the
    recommended approach. Use CVXPYSolver._solve_without_independents() instead.

    Parameters
    ----------
    problem : Problem
        The Problem with relevant matrices and vectors
    solver_convex : str, optional
        Solver to use, by default CLARABEL.
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    tuple
        (fopt, qopt, exitflag, niter, status, sol_time)
    """

    # Retrieve important vector and matrices
    q = problem.q
    E = problem.E
    C = problem.C
    Ci = problem.Ci
    Cit = problem.Cit
    Cb = problem.Cb
    pz = problem.P[:, 2].reshape(-1, 1)
    ph = problem.ph
    free = problem.free
    fixed = problem.fixed
    qmin = problem.qmin
    qmax = problem.qmax
    x = problem.x0
    y = problem.y0
    m = problem.m

    q = cp.Variable(m)
    fobj = matrix_frac(pz[free], -Cit @ cp.diag(q) @ Ci) - x.T @ C.T @ diag(q) @ Cb @ x[fixed] - y.T @ C.T @ diag(q) @ Cb @ y[fixed]  # for q negative
    objective = Minimize(fobj)

    horz = E @ q == ph.flatten()
    pos = q >= qmin.flatten()
    maxq = q <= qmax.flatten()

    constraints = [horz, pos, maxq]

    prob = CVXProblem(objective, constraints)
    prob.solve(solver=solver_convex, verbose=printout)

    # save output
    fopt = prob.value
    qopt = q.value
    status = prob.status
    niter = prob.solver_stats.num_iters
    sol_time = prob.solver_stats.solve_time

    if status not in ["infeasible", "unbounded"]:
        exitflag = 1
    else:
        exitflag = 0

    return fopt, qopt, exitflag, niter, status, sol_time


def call_cvxpy_ind(problem, solver_convex="CLARABEL", printout=False):
    """Legacy helper function - use CVXPYSolver class instead.

    This function is kept for backward compatibility but is no longer the
    recommended approach. Use CVXPYSolver._solve_with_independents() instead.

    Parameters
    ----------
    problem : Problem
        The Problem with relevant matrices and vectors
    solver_convex : str, optional
        Solver to use, by default CLARABEL.
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    tuple
        (fopt, qopt, exitflag, niter, status, sol_time)
    """

    # Retrieve important vector and matrices
    q = problem.q
    Edinv = problem.Edinv
    Ei = problem.Ei
    ph = problem.ph
    C = problem.C
    Ci = problem.Ci
    Cit = problem.Cit
    Cb = problem.Cb
    pz = problem.P[:, 2].reshape(-1, 1)
    free = problem.free
    fixed = problem.fixed
    dep = problem.dep
    ind = problem.ind
    qmin = problem.qmin
    qmax = problem.qmax
    x = problem.x0
    y = problem.y0
    m = problem.m

    q = cp.Variable(m)

    fobj = matrix_frac(pz[free], -Cit @ cp.diag(q) @ Ci) - x.T @ C.T @ diag(q) @ Cb @ x[fixed] - y.T @ C.T @ diag(q) @ Cb @ y[fixed]
    objective = Minimize(fobj)

    horz = q[dep] == Edinv @ (Ei @ q[ind] - ph.flatten())
    pos = q >= qmin.flatten()
    maxq = q <= qmax.flatten()

    constraints = [horz, pos, maxq]

    prob = CVXProblem(objective, constraints)
    prob.solve(solver=solver_convex, verbose=printout)  # open source solver, try MOSEK if CLARABEL fails

    # save output
    fopt = prob.value
    qopt = q.value
    status = prob.status
    niter = prob.solver_stats.num_iters
    sol_time = prob.solver_stats.solve_time

    if status not in ["infeasible", "unbounded"]:
        exitflag = 1
    else:
        exitflag = 0

    return fopt, qopt, exitflag, niter, status, sol_time
