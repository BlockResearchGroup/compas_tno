"""Base solver classes for optimization problems."""

from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING
from typing import Callable
from typing import List
from typing import Optional
from typing import Tuple

import numpy.typing as npt

if TYPE_CHECKING:
    from compas_tno.problems import Problem
    from compas_tno.solvers.result import SolverResult


class Solver(ABC):
    """Abstract base class for optimization solvers.

    This class defines the interface for all solvers in compas_tno.
    Solvers can work with Problem objects and optimization functions.

    Parameters
    ----------
    settings : dict, optional
        Solver-specific settings.

    Examples
    --------
    >>> solver = IPOPTSolver(settings={"printout": True, "max_iter": 500})
    >>> result = solver.solve(problem, objective=fobj, constraints=fconstr)
    >>> if result.success:
    ...     print(f"Optimal value: {result.fopt}")

    """

    def __init__(self, settings: Optional[dict] = None):
        """Initialize the solver.

        Parameters
        ----------
        settings : dict, optional
            Solver-specific settings.

        """
        self.settings = settings or {}

    @abstractmethod
    def solve(self, problem: "Problem", objective: Optional[Callable] = None, constraints: Optional[Callable] = None, **kwargs) -> "SolverResult":
        """Solve the optimization problem.

        Parameters
        ----------
        problem : Problem
            The problem structure with matrices and data.
        objective : callable, optional
            Objective function f(x, problem) -> float.
            If None, uses problem.fobj.
        constraints : callable, optional
            Constraint function g(x, problem) -> array.
            If None, uses problem.fconstr.
        **kwargs
            Additional solver-specific keyword arguments.

        Returns
        -------
        SolverResult
            The optimization result.

        """
        pass


class NonlinearSolver(Solver):
    """Base class for nonlinear optimization solvers.

    This class provides common functionality for nonlinear programming (NLP) solvers
    such as computing starting points and bounds from the problem structure.

    """

    def solve(
        self,
        problem: "Problem",
        objective: Optional[Callable] = None,
        constraints: Optional[Callable] = None,
        gradient: Optional[Callable] = None,
        jacobian: Optional[Callable] = None,
        x0: Optional[npt.NDArray] = None,
        bounds: Optional[List[Tuple[float, float]]] = None,
        callback: Optional[Callable] = None,
        **kwargs,
    ) -> "SolverResult":
        """Solve a nonlinear optimization problem.

        Parameters
        ----------
        problem : Problem
            The problem structure.
        objective : callable, optional
            Objective function. If None, uses problem.fobj.
        constraints : callable, optional
            Constraint function. If None, uses problem.fconstr.
        gradient : callable, optional
            Gradient of objective. If None, uses problem.fgrad.
        jacobian : callable, optional
            Jacobian of constraints. If None, uses problem.fjac.
        x0 : array, optional
            Starting point. If None, computed from problem.
        bounds : list of tuples, optional
            Variable bounds. If None, computed from problem.
        callback : callable, optional
            Callback function called at each iteration.
        **kwargs
            Additional solver-specific options.

        Returns
        -------
        SolverResult
            The optimization result.

        """
        # Get functions from problem if not provided
        obj = objective or problem.fobj
        cons = constraints or problem.fconstr
        grad = gradient or problem.fgrad
        jac = jacobian or problem.fjac

        if obj is None:
            raise ValueError("Objective function must be provided or stored in problem.fobj")

        # Compute x0 and bounds if not provided
        if x0 is None:
            x0 = self._compute_x0(problem)
        if bounds is None:
            bounds = self._compute_bounds(problem)

        # Delegate to specific solver implementation
        return self._solve_nlp(problem=problem, objective=obj, constraints=cons, gradient=grad, jacobian=jac, x0=x0, bounds=bounds, callback=callback, **kwargs)

    @abstractmethod
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
    ) -> "SolverResult":
        """Solver-specific NLP implementation.

        This method must be implemented by concrete solver classes.

        """
        pass

    def _compute_x0(self, problem: "Problem") -> npt.NDArray:
        """Compute starting point from problem.

        Parameters
        ----------
        problem : Problem
            The problem structure.

        Returns
        -------
        array
            Starting point. Uses problem.x0 if available, otherwise
            computes from force densities of independent edges.

        """
        # Use stored x0 if available (from setup)
        if hasattr(problem, "x0") and problem.x0 is not None:
            return problem.x0

        # Fallback: compute from q
        return problem.q[problem.ind].reshape(-1, 1)

    def _compute_bounds(self, problem: "Problem") -> List[Tuple[float, float]]:
        """Compute variable bounds from problem.

        Parameters
        ----------
        problem : Problem
            The problem structure.

        Returns
        -------
        list of tuples
            Bounds for each variable. Uses problem.bounds if available,
            otherwise computes from qmin/qmax of independent edges.

        """
        # Use stored bounds if available (from setup)
        if hasattr(problem, "bounds") and problem.bounds is not None:
            return problem.bounds

        # Fallback: compute from qmin/qmax
        bounds = []
        for i, (qmin, qmax) in enumerate(zip(problem.qmin, problem.qmax)):
            if i in problem.ind:
                bounds.append((qmin.item(), qmax.item()))
        return bounds


class ConvexSolver(Solver):
    """Base class for convex optimization solvers.

    Convex solvers typically have the objective and constraints
    embedded in the problem formulation.

    """

    def solve(self, problem: "Problem", **kwargs) -> "SolverResult":
        """Solve a convex optimization problem.

        Parameters
        ----------
        problem : Problem
            The problem structure.
        **kwargs
            Solver-specific options.

        Returns
        -------
        SolverResult
            The optimization result.

        """
        return self._solve_convex(problem, **kwargs)

    @abstractmethod
    def _solve_convex(self, problem: "Problem", **kwargs) -> "SolverResult":
        """Solver-specific convex implementation."""
        pass
