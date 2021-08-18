from compas_tno.problems import initialize_loadpath


__all__ = ['initialize_loadpath_proxy',
           'run_NLP_proxy'
           ]


def initialize_loadpath_proxy(formdata, problem=None):

    from compas_tno.diagrams import FormDiagram
    form = FormDiagram.from_data(formdata)
    output = initialize_loadpath(form, problem=problem)

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return form.to_data(), output


def run_NLP_proxy(shapedata, formdata, optimiserdata):

    from compas_tno.diagrams import FormDiagram
    from compas_tno.analysis import Analysis
    from compas_tno.optimisers import Optimiser
    from compas_tno.shapes import Shape

    shape = Shape.from_library(shapedata)
    optimiser = Optimiser()
    optimiser.settings = optimiserdata
    form = FormDiagram.from_data(formdata)

    analysis = Analysis.from_elements(shape, form, optimiser)
    analysis.set_up_optimiser()
    analysis.run()

    # TODO: add an option here to initialise the loadpath with SLSQP/IPOPT if matlab/cvx can not be reached.

    return shape.datashape, form.to_data(), optimiser.settings
