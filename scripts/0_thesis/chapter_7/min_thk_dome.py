from compas_tno import analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_bounds_reactions
from compas_tno.problems import initialise_form
from compas.colors import Color

discr = [20, 16]
discr = [12, 12]

path = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_compression_neg.json'

form = FormDiagram.from_json(path)

# form = FormDiagram.create_circular_radial_form(discretisation=discr)
shape = Shape.create_dome(thk=0.20454691527171837)

apply_envelope_from_shape(form, shape)
apply_bounds_reactions(form, shape)

# initialise_form(form)  # compute the independent edges

# # plotter = TNOPlotter(form, shape)
# # plotter.draw_form_independents()
# # plotter.draw_supports(color=Color.red())
# # plotter.show()

# analysis = Analysis.create_minthk_analysis(form, shape)

# analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])

# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()
# analysis.run()

# plotter = TNOPlotter(form, shape)
# plotter.draw_form()
# plotter.draw_cracks()
# plotter.draw_supports()
# plotter.show()

view = Viewer(form, shape=shape)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()
view.show()

objects = view.app.view.objects

from compas_view2.objects import LineObject
from compas_view2.objects import PointObject

data = {}
lines = []
points = []

for obj in objects:
    if isinstance(obj, LineObject):
        print(obj, obj.linewidth, obj.linecolor, obj._data)
        lines.append({'start': [obj._data[0].x, obj._data[0].y, obj._data[0].z],
                      'end': [obj._data[1].x, obj._data[1].y, obj._data[1].z],
                      'color': (255*obj.linecolor).tolist(),
                      'width': obj.linewidth
                      })
    if isinstance(obj, PointObject):
        print(obj, obj.pointsize, obj.pointcolor, obj._data)
        points.append({'pos': [obj._data.x, obj._data.y, obj._data.z],
                       'color': (255*obj.pointcolor).tolist()
                       })

print(lines)
print(points)
