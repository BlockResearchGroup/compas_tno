from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter

form = FormDiagram.create_cross_form()

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.show()
