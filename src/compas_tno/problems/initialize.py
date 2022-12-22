from compas_tno.solvers.solver_MATLAB import run_loadpath_from_form_MATLAB
from compas_tno.solvers.solver_cvxpy import run_loadpath_from_form_CVXPY

from compas_tno.algorithms import form_update_with_parallelisation
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.algorithms import compute_reactions

import importlib.util


def initialize_loadpath(form, problem=None, find_inds=False, solver_convex='CVXPY', printout=False):
    """Built-in function to optimise the loadpath considering diagram fixed projection.
    Note: This function will select the most appropriate solver (CVX or MOSEK)

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to compute the loadpath
    problem : :class:`~compas_tno.problems.Problem`, optional
        The problem classs with the matrices relevant of the problem, by default None
    find_inds : bool, optional
        If independents need to be found before the loadpath computation, by default False
    solver_convex : str, optional
        Solver to compute the convex optimisation, by default CVXPY
    printout : bool, optional
        If prints about the optimisation data should appear in the screen, by default False

    Returns
    -------
    problem : Problem
        The class with the main matrices of the problem
    """

    if solver_convex == 'CVX' or solver_convex == 'MATLAB':
        if not importlib.util.find_spec('matlab'):
            raise ValueError('MATLAB/CVX not configured. Try changing the <solver_convex> attribute.')
        problem = run_loadpath_from_form_MATLAB(form, problem=problem, find_inds=find_inds, printout=printout)
    elif solver_convex == 'CVXPY' or solver_convex == 'MOSEK':
        if not importlib.util.find_spec('cvxpy'):
            raise ValueError('CVXPY/MOSEK not configured. Try changing the <solver_convex> attribute.')
        problem = run_loadpath_from_form_CVXPY(form, problem=problem, find_inds=find_inds, printout=printout)
    else:
        raise ValueError('Could not initilalise loadpath optimisation with {}. Try changing the <solver_convex> attribute.'.format(solver_convex))

    return problem


def initialize_tna(form, plot=False):
    """Initialize the equilibrium in a form diagram with applied loads using TNA interative solver procedure (form and force diagrams are parallel)

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram. Loads and support must already have been assigned
    plot : bool, optional
        Plots of the intermediare and final force diagrams to follow the process, by default False
    """

    form_update_with_parallelisation(form, plot=plot)

    compute_reactions(form)


def initialize_fdm(form):
    """Initialize the equilibrium in a form diagram with applied loads using FDM approach for the q's stored in the form

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram. Loads, supports and force densities must already have been assigned
    """

    equilibrium_fdm(form)

    compute_reactions(form)
