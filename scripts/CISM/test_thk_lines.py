from compas.geometry import Line
from compas.geometry import Point
from compas_view2.app import App

p1 = Point(0, 0, 0)
p2 = Point(2, 2, 0)
p3 = Point(2, 4, 0)

l1 = Line(p1, p2)
l2 = Line(p2, p3)

view = App(show_grid=False)

view.add(l1, linewidth=0.7)
view.add(l2, linewidth=1.3)

view.show()
