from compas_tno.solvers import run_loadpath_from_form_MATLAB


__all__ = ['initialize_loadpath_proxy',
           'test_import_proxy',
           'initialize_loadpath',
           'initialize_tna'
           ]


def test_import_proxy(formdata):

    from compas_tno.diagrams import FormDiagram
    form = FormDiagram.from_data(formdata)

    import ipopt
    print(ipopt)
    print(form)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return


def initialize_loadpath_proxy(formdata, problem=None):

    from compas_tno.diagrams import FormDiagram
    form = FormDiagram.from_data(formdata)
    initialize_loadpath(form, problem=problem)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return


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
