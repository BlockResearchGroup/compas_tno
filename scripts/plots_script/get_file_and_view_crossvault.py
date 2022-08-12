import tkinter as tk
from tkinter import filedialog
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_view2.shapes import Arrow
from compas.geometry import Point
from compas_plotters import Plotter
from compas.geometry import distance_point_point_xy
import math

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename()
print(file_path)

form = FormDiagram.from_json(file_path)

# form = FormDiagram.create_circular_radial_form()
dome = Shape.create_dome(thk=0.50, radius=5.0, discretisation=[50, 50], t=0.0)
dome = Shape.create_crossvault(thk=0.50, xy_span=[[0,10],[0,10]], discretisation=50, t=0.0)

inds = [edge for edge in form.edges_where({'is_ind': True})]
print('Diagram presents {} independents'.format(len(inds)))

# # find closest node to
# xp, yp = 5.0, 5.0
# dist = 0.01
# for key in form.vertices():
#     coords = form.vertex_coordinates(key)
#     if distance_point_point_xy(coords, [xp, yp]) < dist:
#         load_node = key

plotter = TNOPlotter(form)
plotter.draw_form_independents()
plotter.show()

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.show()

plotter = TNOPlotter(form, figsize=(14, 6))
plotter.settings['size.edge.max_thickness'] = 8.0
plotter.draw_form()
plotter.draw_cracks()
plotter.draw_supports()
plotter.draw_force()
plotter.show()

form = FormDiagram.from_json(file_path)

# load_node = 0
# load_node = 0
# load_node = 148  # A2
# load_node = 131  # A2-sym
# load_node = 50  # B2
# load_node = 87  # C2
# load_node = 86  # D2
# load_node = 101  # D2-mod and radial fd 1/4 span
# load_node = 158  # E2

# to view end point
length = 2.0
x, y, z = form.vertex_coordinates(load_node)
z += length + 0.1
arrow = Arrow([x, y, z], [0, 0, -length])

view = Viewer(form, dome)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 30
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 50
view.draw_thrust(absolute_scale=True)
view.draw_shape()
# view.app.add(Point(x, y, z-length), size=10)
view.draw_cracks()
view.app.add(arrow, color=(0, 0, 0))
view.draw_reactions()
view.app.view.camera.rotation.set(math.pi/4, 0, math.pi/4)
view.show()
