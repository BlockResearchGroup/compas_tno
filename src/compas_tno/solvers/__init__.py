from .cvxpy import CVXPYSolver, call_and_output_CVXPY, run_convex_optimisation, run_loadpath_from_form_CVXPY
from .ipopt import IPOPTSolver, run_nlopt_ipopt
from .post_process import apply_solution_to_form, post_process_nlopt
from .scipy import ScipySolver, run_nlopt_scipy
from .solver import ConvexSolver, NonlinearSolver, Solver, SolverResult

__all__ = [
    # New solver classes
    "Solver",
    "SolverResult",
    "NonlinearSolver",
    "ConvexSolver",
    "IPOPTSolver",
    "ScipySolver",
    "CVXPYSolver",
    # New post-processing
    "apply_solution_to_form",
    # Legacy functions (backward compatibility)
    "run_nlopt_scipy",
    "run_nlopt_ipopt",
    "run_convex_optimisation",
    "run_loadpath_from_form_CVXPY",
    "call_and_output_CVXPY",
    "post_process_nlopt",
]
