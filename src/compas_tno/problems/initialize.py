from compas_tno.solvers.solver_MATLAB import run_loadpath_from_form_MATLAB
from compas_tno.solvers.solver_cvxpy import run_loadpath_from_form_CVXPY
import importlib.util


def initialize_loadpath(form, problem=None, find_inds=False, solver_convex='CVX'):
    """Built-in function to optimise the loadpath considering diagram fixed projection.
    Note: This function will select the most appropriate solver (CVX or MOSEK)

    Parameters
    ----------
    form : FormDiagram
        The form diagram to compute the loadpath
    problem : Problem, optional
        The problem classs with the matrices relevant of the problem, by default None
    find_inds : bool, optional
        If independents need to be found before the loadpath computation, by default False
    solver_convex : str, optional
        Solver to compute the convex optimisation, by default CVX

    Returns
    -------
    None
        Optimisation runs in the background
    """

    if solver_convex == 'CVX' or solver_convex == 'MATLAB':
        if not importlib.util.find_spec('matlab'):
            raise ValueError('MATLAB not configured. Try changing the <solver-convex> attribute.')
        run_loadpath_from_form_MATLAB(form, problem=problem, find_inds=find_inds)
    elif solver_convex == 'CVXPY' or solver_convex == 'MOSEK':
        if not importlib.util.find_spec('cvxpy'):
            raise ValueError('MATLAB not configured. Try changing the <solver-convex> attribute.')
        run_loadpath_from_form_CVXPY(form, problem=problem, find_inds=find_inds)


def initialize_tna(form, plot=False):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    from compas_tno.algorithms import form_update_with_parallelisation

    form_update_with_parallelisation(form, plot=plot)
