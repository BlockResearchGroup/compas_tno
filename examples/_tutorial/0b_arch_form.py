from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer

H = 1.0
L = 2.0
b = 0.5
discretisation = 20

arch = Shape.create_arch(H=H, L=L, b=b)

form = FormDiagram.create_arch(H=H, L=L, discretisation=discretisation)

view = Viewer(form, arch)
view.draw_form()
view.draw_shape()
view.show()
