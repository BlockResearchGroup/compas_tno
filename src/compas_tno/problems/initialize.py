from compas_tno.solvers.solver_MATLAB import run_loadpath_from_form_MATLAB


__all__ = ['initialize_loadpath',
           'initialize_tna'
           ]


def initialize_loadpath(form, problem=None):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    output = run_loadpath_from_form_MATLAB(form, problem=problem)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return output


def initialize_tna(form, plot=False):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    from compas_tno.algorithms.graphstatics import form_update_with_parallelisation

    form_update_with_parallelisation(form, plot=plot)

    return
