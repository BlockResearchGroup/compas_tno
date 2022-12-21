from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_view2.shapes import Arrow
from compas.geometry import distance_point_point_xy
from compas.geometry import Point
from compas.colors import Color
import json

discr = [20, 16]
# discr = [12, 12]

folder_mod = '/Users/mricardo/compas_dev/me/max_load/dome/quadspan/modified_diagrams/'

# file = 'dome_added_paths_discr_[16, 20]_max_load_thk_50.0_pct_stw_0.072.json'
# file = 'dome_added_diagonals_discr_[16, 20]_max_load_thk_50.0_pct_stw_0.092.json'
# file = 'dome_new_pole_max_load_thk_50.0_pct_stw_0.099.json'
file = 'dome_new_pole_and_diagonals_max_load_thk_50.0_pct_stw_0.102.json'

# json_file = folder_mod + 'solution_D.json'

path = folder_mod + file

# Normal domes
# path = '/Users/mricardo/compas_dev/me/max_load/dome/quadspan/sensitivity/dome/radial_fd/dome_radial_fd_discr_[16, 20]_max_load_thk_50.0_pct_stw_0.049.json'
# path = '/Users/mricardo/compas_dev/me/max_load/dome/apex/sensitivity/dome/radial_fd/dome_radial_fd_discr_[16, 20]_max_load_thk_50.0_pct_stw_0.144.json'

form = FormDiagram.from_json(path)

inds = 0
edges = 0
supports = 0
max_force = 0
weight = 0

forces = []

for edge in form.edges_where({'_is_edge': True}):
    edges += 1
    f = abs(form.edge_length(*edge) * form.edge_attribute(edge, 'q'))
    forces.append(f)
    if form.edge_attribute(edge, 'is_ind'):
        inds += 1

for key in form.vertices():
    _, _, z = form.vertex_coordinates(key)
    weight += abs(form.vertex_attribute(key, 'pz'))
    if form.vertex_attribute(key, 'is_fixed'):
        supports += 1
        if z > 0.0:
            rx = form.vertex_attribute(key, '_rx')
            ry = form.vertex_attribute(key, '_ry')
            rz = form.vertex_attribute(key, '_rz')
            r =  (rx**2 + ry**2 + rz**2)**(1/2)
            forces.append(r)

max_f = max(forces)

print('Number of edges:', edges)
print('Number of supports:', supports)
print('Number of independents:', inds)
print('Weight:', weight)
print('Max_force', max_f)

shape = Shape.create_dome()
view = Viewer(form, shape=shape)
view.scale_edge_thickness(1.5)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 45  # normally it's 60 deg
# view.settings['camera.rz'] = 0
# view.settings['camera.rx'] = 0
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()

xp, yp = 7.5, 5.0
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    d = distance_point_point_xy([x, y, 0], [xp, yp, 0])
    if d < 1e-2:
        loaded_node = key

length = 2.0
# loaded_node = 0
x, y, z = form.vertex_coordinates(loaded_node)
z += length + 0.1
arrow = Arrow([x, y, z], [0, 0, -length])
view.app.add(arrow, color=(0, 0, 0))

view.show()

plt = TNOPlotter(form)
plt.draw_form(scale_width=True, absolute_scale=True, color=Color.black())
plt.draw_force()
# plt.draw_form_independents()
plt.draw_supports(size=6.0, color=Color.red())
plt.app.add(Point(x, y, z), facecolor=Color.grey(), size=6.0)
plt.show()

print(xxx)

objects = view.app.view.objects

from compas_view2.objects import LineObject
from compas_view2.objects import PointObject

lines = []
points = []

for obj in objects:
    if isinstance(obj, LineObject):
        # print(obj, obj.linewidth, obj.linecolor, obj._data)
        lines.append({'start': [obj._data[0].x, obj._data[0].y, obj._data[0].z],
                      'end': [obj._data[1].x, obj._data[1].y, obj._data[1].z],
                      'color': (255*obj.linecolor).tolist(),
                      'width': obj.linewidth
                      })
    if isinstance(obj, PointObject):
        # print(obj, obj.pointsize, obj.pointcolor, obj._data)
        points.append({'pos': [obj._data.x, obj._data.y, obj._data.z],
                       'color': (255*obj.pointcolor).tolist()
                       })

# print(lines)
# print(points)

data = {}
data['lines'] = lines
data['points'] = points

with open(json_file, "w") as outfile:
    json.dump(data, outfile)

print('Saved in:', json_file)
