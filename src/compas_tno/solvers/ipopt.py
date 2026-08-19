import time
from typing import TYPE_CHECKING
from typing import Callable
from typing import List
from typing import Optional
from typing import Tuple

import numpy.typing as npt

from compas_tno.solvers.result import SolverResult
from compas_tno.solvers.solver import NonlinearSolver

try:
    import cyipopt

    HAS_IPOPT = True
except ImportError:
    HAS_IPOPT = False

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.problems import Problem


class Wrapper_ipopt:
    """Wrapper to send to IPOPT using analytical derivatives"""

    def __init__(self):
        self.fobj = None
        self.fconstr = None
        self.args = None
        self.fgrad = None
        self.bounds = None
        self.callback = None
        self.x0 = None
        self.eps = 1e-8
        self.fjac = None
        pass

    def objective(self, x):
        """The callback for calculating the objective

        Parameters
        ----------
        x : array
            The variabless passed to IPOPT

        Returns
        -------
        fopt : float
            The objective function value at x.
        """

        if self.callback:
            self.callback(x, *self.args)
        return self.fobj(x, *self.args)

    def gradient(self, x):
        """The callback for calculating the gradient

        Parameters
        ----------
        x : array
            The variabless passed to IPOPT

        Returns
        -------
        grad : array
            The gradient of the objective function at x.
        """

        return self.fgrad(x, *self.args)

    def constraints(self, x):
        """The callback for calculating the constraints

        Parameters
        ----------
        x : array
            The variabless passed to IPOPT

        Returns
        -------
        constr : array
            The constraints of the objective function at x.
        """

        return self.fconstr(x, *self.args).reshape(-1, 1)

    def jacobian(self, x):
        """The callback for calculating the jacobian

        Parameters
        ----------
        x : array
            The variabless passed to IPOPT

        Returns
        -------
        jac : array
            The gradient of the jacobian matrix at x.
        """

        return self.fjac(x, *self.args).flatten()


class IPOPTSolver(NonlinearSolver):
    """IPOPT solver for nonlinear optimization problems.

    This solver uses the Interior Point OPTimizer (IPOPT) for solving
    large-scale nonlinear programming problems.

    Parameters
    ----------
    settings : dict, optional
        Solver settings. Common options:
        - 'printout' : bool - Print solver output (default: False)
        - 'max_iter' : int - Maximum iterations (default: 3000)
        - 'derivative_test' : bool - Test derivatives (default: False)
        - 'nlp_scaling_method' : str - Scaling method (default: None)

    Examples
    --------
    Minimal usage - everything from problem:

    >>> problem.fobj = my_objective
    >>> problem.fconstr = my_constraints
    >>> problem.fgrad = my_gradient
    >>> problem.fjac = my_jacobian
    >>> problem.x0 = initial_guess
    >>> problem.bounds = variable_bounds
    >>> solver = IPOPTSolver(settings={"printout": True})
    >>> result = solver.solve(problem)  # All data from problem!

    Override specific functions:

    >>> result = solver.solve(problem, objective=custom_obj)  # Custom objective only

    """

    def __init__(self, settings: Optional[dict] = None):
        """Initialize IPOPT solver."""
        if not HAS_IPOPT:
            raise ImportError("IPOPT is not installed. Install with: pip install cyipopt")
        super().__init__(settings)

    def _solve_nlp(
        self,
        problem: "Problem",
        objective: Callable,
        constraints: Optional[Callable],
        gradient: Optional[Callable],
        jacobian: Optional[Callable],
        x0: npt.NDArray,
        bounds: List[Tuple[float, float]],
        callback: Optional[Callable],
        **kwargs,
    ) -> SolverResult:
        """Solve using IPOPT.

        Parameters
        ----------
        problem : Problem
            The problem structure.
        objective : callable
            Objective function f(x, problem) -> float.
        constraints : callable, optional
            Constraint function g(x, problem) -> array.
        gradient : callable, optional
            Gradient function grad_f(x, problem) -> array.
        jacobian : callable, optional
            Jacobian function jac_g(x, problem) -> array.
        x0 : array
            Starting point.
        bounds : list of tuples
            Variable bounds.
        callback : callable, optional
            Callback function.
        **kwargs
            Additional options.

        Returns
        -------
        SolverResult
            The optimization result.

        """
        printout = self.settings.get("printout", False)
        args = [problem]

        # Prepare bounds
        lower = [lw for lw, _ in bounds]
        upper = [up for _, up in bounds]

        # Create wrapper for IPOPT
        problem_obj = Wrapper_ipopt()
        problem_obj.fobj = objective
        problem_obj.fconstr = constraints
        problem_obj.fjac = jacobian
        problem_obj.args = args
        problem_obj.fgrad = gradient
        problem_obj.bounds = bounds
        problem_obj.callback = callback
        problem_obj.x0 = x0

        # Compute initial constraints
        if constraints:
            g0 = constraints(x0, *args)
            cu = [10e10] * len(g0)
            cl = [0.0] * len(g0)
            m = len(g0)
        else:
            cu = []
            cl = []
            m = 0

        # Create IPOPT problem
        nlp = cyipopt.Problem(n=len(x0), m=m, problem_obj=problem_obj, lb=lower, ub=upper, cl=cl, cu=cu)

        # Set options
        nlp = self._set_options(nlp)

        # Solve
        start_time = time.time()
        xopt, info = nlp.solve(x0.flatten())
        elapsed_time = time.time() - start_time

        # Extract results
        f_opt = info["obj_val"]
        exit_flag = info["status"]

        # IPOPT: status=0 or 1 is success
        if exit_flag == 1 or exit_flag == 0:
            success = True
            exit_flag = 0
            message = "Solved successfully with IPOPT"
        else:
            success = False
            exit_flag = 1
            message = "Error: Did not find convergence (IPOPT)"

        if printout:
            print(str(info["status_msg"]))
            print(f"Solving Time: {elapsed_time:.1f} sec")

        return SolverResult(
            xopt=xopt.reshape(-1, 1),
            fopt=float(f_opt),
            success=success,
            message=message,
            niter=None,  # IPOPT doesn't easily expose iteration count
            time=elapsed_time,
            exitflag=exit_flag,
        )

    def _set_options(self, nlp) -> "cyipopt.Problem":
        """Set IPOPT solver options.

        Parameters
        ----------
        nlp : cyipopt.Problem
            The IPOPT problem object.

        Returns
        -------
        cyipopt.Problem
            The problem with options set.

        """
        if not self.settings.get("printout", False):
            nlp.add_option("print_level", 0)

        if self.settings.get("derivative_test", False):
            nlp.add_option("derivative_test", "first-order")

        if self.settings.get("max_iter"):
            nlp.add_option("max_iter", self.settings["max_iter"])

        scaling = self.settings.get("nlp_scaling_method")
        if scaling:
            if self.settings.get("printout"):
                print(f"Applied solver scaling: {scaling}")
            nlp.add_option("nlp_scaling_method", scaling)

        return nlp


def run_nlopt_ipopt(analysis: "Analysis"):
    """Run nonlinear optimisation problem with IPOPT.

    This function is kept for backward compatibility. New code should use IPOPTSolver class.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    result : Result
        The Result object with the key solution information.

    Notes
    -----
    This function wraps the new IPOPTSolver class for backward compatibility.
    For new code, prefer using::

        solver = IPOPTSolver(settings=optimiser.settings)
        result = solver.solve(problem, objective=fobj, constraints=fconstr, ...)

    """
    optimiser = analysis.optimiser
    problem = optimiser.problem

    # Create solver using optimiser settings
    solver = IPOPTSolver(settings=optimiser.settings)

    # Solve using the new class - functions and variables come from problem
    result = solver.solve(problem)

    # # Store results in optimiser
    # optimiser.exitflag = result.exitflag
    # optimiser.time = result.time
    # optimiser.fopt = result.fopt
    # optimiser.xopt = result.xopt
    # optimiser.niter = result.niter
    # optimiser.message = result.message

    return result


def _nlp_options(nlp, optimiser: "Optimiser"):
    """Set NLP options for IPOPT

    Parameters
    ----------
    nlp : obj
        IPOPT object
    optimiser : Optimiser
        The Optimiser object

    Returns
    -------
    [type]
        [description]
    """

    # Link to instructions: https://coin-or.github.io/Ipopt/OPTIONS.html

    # nlp.add_option(b'hessian_approximation', b'limited-memory')
    # nlp.add_option('tol', 1e-4)              # Default 1e-8
    # nlp.add_option('max_iter', 500)        # Default 3000

    # nlp.add_option('dual_inf_tol', 100.0)  # Default 1.0
    # nlp.add_option('constr_viol_tol', 0.1) # Default 1e-4
    # nlp.add_option('compl_inf_tol', 0.1)   # Default 1e-4

    # nlp.add_option('acceptable_iter', 10)
    # nlp.add_option('acceptable_tol', 1e-4)  # Default 1e-6
    # nlp.add_option('acceptable_constr_viol_tol', 1e-4)  # Default 1e-2
    # nlp.add_option('acceptable_dual_inf_tol', 10e10)  # Default 10e10
    # nlp.add_option('acceptable_compl_inf_tol', 1e-2)  # Default 1e-2
    # nlp.add_option('max_iter', 500)

    scaling = optimiser.settings.get("nlp_scaling_method", None)

    if not optimiser.settings["printout"]:
        nlp.add_option("print_level", 0)

    if optimiser.settings.get("derivative_test", None):
        nlp.add_option("derivative_test", "first-order")
        # nlp.add_option('derivative_test_perturbation') #, 10e-8
        # nlp.add_option('derivative_test_print_all', 'yes')

    if optimiser.settings.get("max_iter", None):
        nlp.add_option("max_iter", optimiser.settings["max_iter"])

    if scaling:
        print("Applied solver scaling:", scaling)
        nlp.add_option("nlp_scaling_method", scaling)

    return nlp
