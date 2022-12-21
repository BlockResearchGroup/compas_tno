from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer

H = 1.0
L = 2.0
b = 0.5

arch = Shape.create_arch(H=H, L=L, b=b)

view = Viewer(shape=arch, show_grid=True)
view.draw_shape()
view.show()
