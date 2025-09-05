import time
from typing import TYPE_CHECKING

try:
    import cyipopt

    HAS_IPOPT = True
except ImportError:
    HAS_IPOPT = False

if TYPE_CHECKING:
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser


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


def run_nlopt_ipopt(analysis: "Analysis"):
    """Run nonlinear optimisation problem with IPOPT

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis
        Analysis object optimised.

    """

    if not HAS_IPOPT:
        raise ImportError("IPOPT is not installed. Please install it using `pip install cyipopt`")

    optimiser = analysis.optimiser
    printout = optimiser.settings.get("printout", False)
    callback = optimiser.callback

    bounds = optimiser.bounds
    x0 = optimiser.x0
    g0 = optimiser.g0
    args = [optimiser.problem]

    lower = [lw[0] for lw in bounds]
    upper = [up[1] for up in bounds]

    problem_obj = Wrapper_ipopt()
    problem_obj.fobj = optimiser.fobj
    problem_obj.fconstr = optimiser.fconstr
    problem_obj.fjac = optimiser.fjac
    problem_obj.args = args
    problem_obj.fgrad = optimiser.fgrad
    problem_obj.bounds = bounds
    problem_obj.callback = callback
    problem_obj.x0 = x0

    if printout:
        g0 = optimiser.fconstr(x0, *args)

    cu = [10e10] * len(g0)
    cl = [0.0] * len(g0)

    nlp = cyipopt.Problem(n=len(x0), m=len(g0), problem_obj=problem_obj, lb=lower, ub=upper, cl=cl, cu=cu)

    # Set Options and Time
    nlp = _nlp_options(nlp, optimiser)
    start_time = time.time()

    # Solve
    xopt, info = nlp.solve(x0)
    fopt = info["obj_val"]
    exitflag = info["status"]
    if exitflag == 1 or exitflag == 0:  # IPOPT consider solved = 1. Solved in tolerances 0 -> TNO solved = 0
        exitflag = 0
        msg = "Solved successfully with IPOPT"
    else:
        exitflag = 1
        msg = "Error: Did not find convergence (IPOPT)"
    if printout:
        print(str(info["status_msg"]))

    elapsed_time = time.time() - start_time

    if printout:
        print("Solving Time: {0:.1f} sec".format(elapsed_time))

    optimiser.exitflag = exitflag
    optimiser.time = elapsed_time
    optimiser.fopt = fopt
    optimiser.xopt = xopt
    optimiser.niter = None  # Did not find a way to display number of iterations
    # optimiser.message = str(info['status_msg'])
    optimiser.message = msg  # temporary, should be equal to the line above
    optimiser.nlp = nlp

    return analysis


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
        print("Applied sollver scaling:", scaling)
        nlp.add_option("nlp_scaling_method", scaling)

    return nlp
