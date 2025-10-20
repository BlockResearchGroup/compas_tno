"""Base solver classes for optimization problems."""

from dataclasses import dataclass
from typing import Optional

import numpy.typing as npt


@dataclass
class SolverResult:
    """Result from an optimization solver.

    Attributes
    ----------
    xopt : array
        Optimal solution vector.
    fopt : float
        Optimal objective function value.
    success : bool
        Whether the solver succeeded.
    message : str
        Status message from the solver.
    niter : int
        Number of iterations.
    time : float
        Time taken to solve (seconds).
    exitflag : int
        Exit flag (0 = success, non-zero = failure/warning).

    """

    xopt: npt.NDArray
    fopt: float
    success: bool
    message: str
    niter: Optional[int] = None
    time: float = 0.0
    exitflag: int = 0
