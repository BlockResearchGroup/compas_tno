from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_view2.shapes import Arrow
from compas.geometry import distance_point_point_xy
from compas.geometry import Point
from compas.colors import Color
import json

type_diag = 'radial_fd'
discretisation = [16, 20]
# discr = [12, 12]

folder = '/Users/mricardo/compas_dev/me/hor-loads/dome/sensitivity/'
title = 'dome_hor_load_{}_discr_{}.json'.format(type_diag, discretisation)
path = folder + title

form = FormDiagram.from_json(path)

lamdh = form.attributes['lambdh']
print('Lambd:', lamdh)

shape = Shape.create_dome()
view = Viewer(form, shape=shape)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 45  # normally it's 60 deg
view.scale_edge_thickness(2.0)
# view.settings['camera.rz'] = 0
# view.settings['camera.rx'] = 0
view.draw_form(cull_negative=True, absolute_scale=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()

view.show()

plt = TNOPlotter(form)
plt.draw_form(scale_width=False, color=Color.black())
plt.draw_supports(size=6.0, color=Color.red())
plt.show()

plt = TNOPlotter(form)
plt.show_solution()
plt.show()

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

json_file = '/Users/mricardo/compas_dev/me/hor-loads/dome/sensitivity/dome_hor_load_radial_fd_discr_[16, 20]_raw.json'

with open(json_file, "w") as outfile:
    json.dump(data, outfile)

print('Saved in:', json_file)
