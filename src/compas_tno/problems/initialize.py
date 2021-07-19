from compas_tno.solvers import run_loadpath_from_form_MATLAB


__all__ = ['initialize_loadpath',
           'initialize_tna'
           ]


def initialize_loadpath(form, problem=None):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    run_loadpath_from_form_MATLAB(form, problem=problem)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return


def initialize_tna(form, problem=None):
    """ Built-in function to optimise the loadpath considering diagram fixed projection
    """

    print('WIP')

    return
