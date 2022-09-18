import compas_tno
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities import apply_envelope_from_shape

discretisation = 14
thk = 0.5
span = 10.0
spr_angle = 30.0
alpha = 1/math.cos(math.radians(spr_angle))
L = span * alpha
Ldiff = L - span
xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]
shape = Shape.create_crossvault(xy_span=xyspan_shape, thk=thk)

path = compas_tno.open_dialog()\
form = FormDiagram.from_json(path)

plot = TNOPlotter(form)
plot.draw_form_independents()
plot.draw_supports()
plot.show()

base = FormDiagram.create_cross_form(discretisation=50)
apply_envelope_from_shape(base, shape)
shape_fixed = Shape.from_formdiagram_and_attributes(base)

view = Viewer(form, shape_fixed)
view.scale_edge_thickness(10.0)
view.draw_thrust()
view.draw_shape()
view.draw_reactions()
view.draw_cracks()
view.show()
