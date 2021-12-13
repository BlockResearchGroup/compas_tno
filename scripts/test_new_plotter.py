from compas_plotters.artists import segmentartist
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters2 import FormPlotter
import compas_tno

from compas.geometry import Vector
from compas.geometry import Rotation
import math

json = compas_tno.get('form-test.json')

form = FormDiagram.from_json(json)

# plotter = FormPlotter(form)
# plotter.draw_form()
# plotter.draw_cracks()
# plotter.draw_reactions()
# plotter.show()

# axis = Vector(1.0, 0, 0)
# form_rotated = form.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))

plotter = FormPlotter(form)
plotter.draw_form_xz()
plotter.draw_cracks()
plotter.show()
