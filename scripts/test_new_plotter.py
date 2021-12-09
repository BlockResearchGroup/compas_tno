from compas_plotters.artists import segmentartist
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters2 import FormPlotter
import compas_tno

json = compas_tno.get('form-test.json')

form = FormDiagram.from_json(json)

plotter = FormPlotter(form)
plotter.draw_form()
plotter.draw_cracks()
plotter.draw_reactions()
plotter.show()
