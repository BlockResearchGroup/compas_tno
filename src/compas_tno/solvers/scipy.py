import time
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable
from typing import List
from typing import Optional
from typing import Tuple

import numpy as np
import numpy.typing as npt
from scipy.optimize import fmin_slsqp

from compas_tno.solvers.solver import NonlinearSolver
from compas_tno.solvers.solver import SolverResult

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis
    from compas_tno.problems import Problem


class ScipySolver(NonlinearSolver):
    """SciPy SLSQP solver for nonlinear optimization problems.

    This solver uses Sequential Least SQuares Programming (SLSQP)
    from scipy.optimize for constrained nonlinear optimization.

    Parameters
    ----------
    settings : dict, optional
        Solver settings. Common options:
        - 'printout' : bool - Print solver output (default: True)
        - 'max_iter' : int - Maximum iterations (default: 500)
        - 'gradient' : bool - Use analytical gradient (default: False)
        - 'jacobian' : bool - Use analytical jacobian (default: False)

    Examples
    --------
    Minimal usage - everything from problem:

    >>> problem.fobj = my_objective
    >>> problem.fconstr = my_constraints
    >>> problem.fgrad = my_gradient
    >>> problem.fjac = my_jacobian
    >>> problem.x0 = initial_guess
    >>> problem.bounds = variable_bounds
    >>> solver = ScipySolver(settings={"printout": True})
    >>> result = solver.solve(problem)  # All data from problem!

    Override specific functions:

    >>> result = solver.solve(problem, gradient=custom_grad)  # Custom gradient only

    """

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
        """Solve using SciPy SLSQP.

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
        printout = self.settings.get("printout", True)
        grad_choice = self.settings.get("gradient", False)
        jac_choice = self.settings.get("jacobian", False)
        max_iter = self.settings.get("max_iter", 500)

        # Disable gradient/jacobian if settings say so
        if not grad_choice:
            gradient = None
        if not jac_choice:
            jacobian = None

        # Flatten gradient for SLSQP compatibility
        if gradient:
            gradient = self._slsqp_grad_flatten(gradient)

        args = [problem]

        start_time = time.time()

        # Call SLSQP
        fopt, xopt, exitflag, niter, message = self._slsqp(
            fn=objective, qid0=x0, bounds=bounds, fprime=gradient, fprime_ieqcons=jacobian, printout=printout, fieq=constraints, args=args, iter=max_iter, callback=callback
        )

        elapsed_time = time.time() - start_time

        if printout:
            print(f"Solving Time: {elapsed_time:.1f} sec")

        return SolverResult(xopt=xopt.reshape(-1, 1), fopt=float(fopt), success=(exitflag == 0), message=message, niter=niter, time=elapsed_time, exitflag=exitflag)

    def _slsqp(
        self,
        fn: Callable[[np.ndarray, Any], float],
        qid0: np.ndarray,
        bounds: List[Tuple[float, float]],
        fprime: Optional[Callable[[np.ndarray, Any], np.ndarray]],
        fprime_ieqcons: Optional[Callable[[np.ndarray, Any], np.ndarray]],
        printout: bool,
        fieq: Optional[Callable[[np.ndarray, Any], np.ndarray]],
        args: List[Any],
        iter: int,
        callback: Optional[Callable[[np.ndarray], None]],
    ) -> Tuple[float, np.ndarray, int, int, str]:
        """Run SLSQP optimization.

        Parameters
        ----------
        fn : callable
            Objective function.
        qid0 : array
            Starting point.
        bounds : list of tuples
            Variable bounds.
        fprime : callable, optional
            Gradient of objective.
        fprime_ieqcons : callable, optional
            Jacobian of constraints.
        printout : bool
            Whether to print output.
        fieq : callable, optional
            Inequality constraint function.
        args : list
            Additional arguments.
        iter : int
            Maximum iterations.
        callback : callable, optional
            Callback function.

        Returns
        -------
        tuple
            (fopt, xopt, exitflag, niter, message)

        """
        pout = 2 if printout else 0
        opt = fmin_slsqp(fn, qid0, args=args, disp=pout, fprime=fprime, f_ieqcons=fieq, fprime_ieqcons=fprime_ieqcons, bounds=bounds, full_output=1, iter=iter, callback=callback)
        # Return order: fopt, xopt, exitflag, niter, message
        return opt[1], opt[0], opt[3], opt[2], opt[4]

    def _slsqp_grad_flatten(self, fgrad: Callable) -> Callable:
        """Flatten the gradient function for SLSQP compatibility.

        Parameters
        ----------
        fgrad : callable
            Gradient function.

        Returns
        -------
        callable
            Flattened gradient function.

        """

        def _fgrad_flatten(x: np.ndarray, *args: Any) -> np.ndarray:
            grad = fgrad(x, *args)
            return grad.flatten()

        return _fgrad_flatten


def run_nlopt_scipy(analysis: "Analysis") -> "Analysis":
    """Run nonlinear optimisation problem with SciPy.

    This function is kept for backward compatibility. New code should use ScipySolver class.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : Analysis
        Analysis object optimised.

    Notes
    -----
    This function wraps the new ScipySolver class for backward compatibility.
    For new code, prefer using::

        solver = ScipySolver(settings=optimiser.settings)
        result = solver.solve(problem, objective=fobj, constraints=fconstr, ...)

    """
    optimiser = analysis.optimiser
    problem = optimiser.problem

    # Create solver using optimiser settings
    solver = ScipySolver(settings=optimiser.settings)

    # Solve using the new class - functions and variables come from problem
    result = solver.solve(problem)

    # Store results in optimiser
    optimiser.exitflag = result.exitflag
    optimiser.time = result.time
    optimiser.fopt = result.fopt
    optimiser.xopt = result.xopt
    optimiser.niter = result.niter
    optimiser.message = result.message

    return analysis
