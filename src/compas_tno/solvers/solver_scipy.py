import time
from typing import TYPE_CHECKING
from typing import Any, Callable, List, Optional, Tuple, Dict

import numpy as np
from scipy.optimize import fmin_slsqp
from scipy.optimize import shgo

from .post_process import post_process_general

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.problems import Problem

def run_optimisation_scipy(analysis: "Analysis") -> "Analysis":
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
    solver: str = optimiser.settings["solver"]
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

    if solver == "slsqp" or solver == "SLSQP":
        fopt: float
        xopt: np.ndarray
        exitflag: int
        niter: int
        message: str
        fopt, xopt, exitflag, niter, message = _slsqp(
            fobj, x0, bounds, slsqp_grad_flatten(fgrad) if fgrad else None, fjac,
            printout, fconstr, args, max_iter, callback
        )
    elif solver == "shgo":
        dict_constr: List[Dict[str, Any]] = []
        for i in range(len(fconstr(x0, *args))):
            args_constr: List[Any] = list(args)
            args_constr.append(i)
            args_constr.append(fconstr)
            dict_: Dict[str, Any] = {
                "type": "ineq",
                "fun": _shgo_constraint_wrapper,
                "args": args_constr,
            }
            dict_constr.append(dict_)
        args_constr[len(args_constr) - 2] = 0  # type: ignore
        result = _shgo(fobj, bounds, True, dict_constr, args)
        fopt = result["fun"]
        xopt = result["x"]
        sucess = result["success"]
        message = result["message"]
        if sucess is True:
            exitflag = 0
        else:
            exitflag = 1
            print(message)
        niter = None  # Not available for shgo
    else:
        raise ValueError(f"Unknown solver: {solver}")

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

    post_process_general(analysis)

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
    pout = 2 if printout else 0
    opt = fmin_slsqp(
        fn, qid0, args=args, disp=pout, fprime=fprime, f_ieqcons=fieq,
        fprime_ieqcons=fprime_ieqcons, bounds=bounds, full_output=1, iter=iter, callback=callback
    )
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

def _shgo(
    fn: Callable[[np.ndarray, Any], float],
    bounds: List[Tuple[float, float]],
    printout: bool,
    dict_constr: List[Dict[str, Any]],
    args: List[Any]
) -> Dict[str, Any]:
    res = shgo(
        fn, bounds, args=args, constraints=dict_constr,
        n=3, iters=3, options={"disp": True}
    )
    return res

def _shgo_constraint_wrapper(
    x: np.ndarray,
    *args_constr: Any
) -> float:
    length = len(args_constr)
    args = args_constr[: length - 2]
    i = args_constr[length - 2]
    fconstr = args_constr[length - 1]
    gi = fconstr(x, *args)[i]
    return gi

def _cobyla(
    fn: Callable,
    qid0: np.ndarray,
    bounds: List[Tuple[float, float]],
    printout: bool,
    fieq: Optional[Callable],
    args: List[Any]
) -> Optional[Any]:
    # pout = 2 if printout else 0
    # opt  = fmin_cobyla(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=500, f_ieqcons=fieq)

    # W.I.P

    return None

# def _diff_evo(fn, bounds, population, generations, printout, plot, frange, args):
#     return devo_numpy(fn=fn, bounds=bounds, population=population, generations=generations, printout=printout, plot=plot, frange=frange, args=args)

# def _ga(fn, fit_type, num_var, boundaries, num_gen, num_pop, args):
#     return ga(fit_function=fn, fit_type=fit_type, num_var=num_var, boundaries=boundaries, num_gen=num_gen, num_pop=num_pop, fargs=args)
