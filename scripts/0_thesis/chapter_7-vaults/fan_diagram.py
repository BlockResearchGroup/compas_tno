from compas_tno.diagrams import FormDiagram
from compas_tno.problems import initialise_form
from compas_tno.plotters import TNOPlotter

# form = FormDiagram.create_fan_form(discretisation=14)
# initialise_form(form)

# plot = TNOPlotter(form)
# plot.draw_form_independents()
# plot.draw_supports()
# plot.show()


for d in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    form = FormDiagram.create_parametric_form(lambd=d, discretisation=14)
    # initialise_form(form)

    plot = TNOPlotter(form)
    plot.draw_form_independents()
    plot.draw_supports()
    plot.show()
