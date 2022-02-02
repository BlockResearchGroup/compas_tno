from compas_tno.solvers.solver_MATLAB import run_loadpath_from_form_MATLAB
from compas_tno.solvers.solver_cvxpy import run_loadpath_from_form_CVXPY


def initialize_loadpath(form, problem=None, find_inds=False):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    output = run_loadpath_from_form_MATLAB(form, problem=problem, find_inds=find_inds)

    return output


def initialize_loadpath_no_matlab(form, problem=None, find_inds=False):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    output = run_loadpath_from_form_CVXPY(form, problem=problem, find_inds=find_inds)

    return output


def initialize_tna(form, plot=False):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    from compas_tno.algorithms.graphstatics import form_update_with_parallelisation

    form_update_with_parallelisation(form, plot=plot)

    return
