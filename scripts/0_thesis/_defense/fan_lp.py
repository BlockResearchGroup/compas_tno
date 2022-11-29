from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas.colors import Color

form = FormDiagram.create_fan_form(discretisation=14)
form = FormDiagram.create_cross_form(discretisation=14)
form = FormDiagram.create_parametric_form(discretisation=14, lambd=0.5)
form = FormDiagram.create_delta_form(discretisation=14, delta=0.5)

plotter = TNOPlotter(form)

plotter.draw_form(scale_width=False)
plotter.draw_supports(color=Color.red(), size=8)
plotter.show()
