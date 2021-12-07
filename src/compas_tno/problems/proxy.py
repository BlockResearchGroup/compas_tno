from compas_tno.problems import initialize_loadpath


def initialize_loadpath_proxy(formdata, problem=None):
    """ Initialise loadpath using the proxy."""

    from compas_tno.diagrams import FormDiagram

    form = FormDiagram.from_data(formdata)
    output = initialize_loadpath(form, problem=problem)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return form.to_data(), output


def run_NLP_proxy(shapedata, formdata, optimiserdata):
    """ Run nonlinear multiobjective optimisations using the proxy."""

    from compas_tno.diagrams import FormDiagram
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.shapes import Shape

    shape = Shape.from_data(shapedata)
    form = FormDiagram.from_data(formdata)
    optimiser = Optimiser.from_data(optimiserdata)

    analysis = Analysis.from_elements(shape, form, optimiser)
    analysis.set_up_optimiser()
    analysis.run()

    return shape.to_data(), form.to_data(), optimiser.to_data()
