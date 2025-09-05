import time
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable
from typing import List
from typing import Optional
from typing import Tuple

import numpy as np
from scipy.optimize import fmin_slsqp

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.problems import Problem


def run_nlopt_scipy(analysis: "Analysis") -> "Analysis":
    """Run nonlinear optimisation problem with SciPy.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : Analysis
        Analysis object optimised.
    """

    optimiser: "Optimiser" = analysis.optimiser
    fobj: Callable = optimiser.fobj
    fconstr: Callable = optimiser.fconstr
    fgrad: Optional[Callable] = optimiser.fgrad
    fjac: Optional[Callable] = optimiser.fjac
    args: List["Problem"] = [optimiser.problem]
    bounds: List[Tuple[float, float]] = optimiser.bounds
    x0: np.ndarray = optimiser.x0
    printout: bool = optimiser.settings.get("printout", True)
    grad_choice: bool = optimiser.settings.get("gradient", False)
    jac_choice: bool = optimiser.settings.get("jacobian", False)
    max_iter: int = optimiser.settings.get("max_iter", 500)
    callback: Optional[Callable] = optimiser.callback

    if grad_choice is False:
        fgrad = None
    if jac_choice is False:
        fjac = None

    start_time = time.time()

    fopt: float
    xopt: np.ndarray
    exitflag: int
    niter: int
    message: str
    fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, slsqp_grad_flatten(fgrad) if fgrad else None, fjac, printout, fconstr, args, max_iter, callback)

    elapsed_time: float = time.time() - start_time
    if printout:
        print("Solving Time: {0:.1f} sec".format(elapsed_time))

    # Store output info in optimiser
    optimiser.exitflag = exitflag
    optimiser.time = elapsed_time
    optimiser.fopt = float(fopt)
    optimiser.xopt = xopt
    optimiser.niter = niter
    optimiser.message = message

    return analysis


def _slsqp(
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
    """Run SLSQP optimisation."""

    pout = 2 if printout else 0
    opt = fmin_slsqp(fn, qid0, args=args, disp=pout, fprime=fprime, f_ieqcons=fieq, fprime_ieqcons=fprime_ieqcons, bounds=bounds, full_output=1, iter=iter, callback=callback)
    # Return order: fopt, xopt, exitflag, niter, message
    return opt[1], opt[0], opt[3], opt[2], opt[4]


def slsqp_grad_flatten(fgrad: Optional[Callable]) -> Optional[Callable]:
    """Flatten the gradient function to be compatible with SLSQP.

    Parameters
    ----------
    fgrad : Callable
        Gradient function.

    Returns
    -------
    Callable
        Flattened gradient function.
    """
    if fgrad is None:
        return None

    def _fgrad_flatten(x: np.ndarray, *args: Any) -> np.ndarray:
        grad = fgrad(x, *args)
        return grad.flatten()

    return _fgrad_flatten
