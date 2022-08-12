from arrow import Arrow
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
from compas_tno.problems import initialise_form
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.utilities import apply_envelope_from_shape
from compas.colors import Color

# Figures of the initial solutions
path = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-general1-Ecomp-linear-10-thk-0.5-corner-diagonal.json'

# form = FormDiagram.from_json('//Users/mricardo/compas_dev/compas_tno/data/CISM/form-max-16.json')
# form = FormDiagram.from_json('//Users/mricardo/compas_dev/compas_tno/data/CISM/form-min-16.json')
# form = FormDiagram.from_json('//Users/mricardo/compas_dev/compas_tno/data/CISM/form-t-16.json')
# form = FormDiagram.from_json('//Users/mricardo/compas_dev/compas_tno/data/CISM/form-direct_path-max_load-16-thk-0.5.json')
# form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-temp.json')
form = FormDiagram.from_json(path)

swt = form.lumped_swt()
thrust = form.thrust()
print('T/W:', round(thrust/swt, 2))

vault = Shape.from_formdiagram_and_attributes(form)

# Plot form force
plotter = TNOPlotter(form, figsize=(14, 6))
plotter.settings['size.edge.max_thickness'] = 8.0
plotter.settings['color.edges.form'] = Color.black()
plotter.settings['color.vertex.supports'] = Color.red()
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.draw_force()
plotter.show()

# from compas_plotters import Plotter
# plotter = Plotter()
# artist = plotter.add(form)
# artist.draw_vertexlabels()
# plotter.show()

# Classic plot
plotter = TNOPlotter(form)
plotter.show_solution()

viewer = Viewer(form, shape=vault, show_grid=False)
viewer.settings['camera.target'] = [5, 5, 0]
viewer.settings['camera.distance'] = 35
viewer.settings['scale.reactions'] = 0.005
viewer.settings['scale.reactions'] = 0.005/3
viewer.settings['scale.reactions'] = 0.005/6
# viewer.settings['scale.reactions'] = 0.005*2
viewer.settings['opacity.shapes'] =  0.3
viewer.draw_thrust()

# from compas_view2.objects import Arrow
# loaded_node = 143
# length = 2.0
# x, y, z = form.vertex_coordinates(loaded_node)
# z += length + 0.1
# arrow = Arrow([x, y, z], [0, 0, -length])
# viewer.app.add(arrow, color=(0, 0, 0))

viewer.draw_cracks()
viewer.draw_shape()
viewer.draw_reactions()
viewer.show()

# Classic plot
plotter = TNOPlotter(form)
plotter.show_solution()

viewer = Viewer(form, shape=vault, show_grid=False)
viewer.settings['camera.target'] = [5, 5, 0]
viewer.settings['camera.distance'] = 35
viewer.settings['scale.reactions'] = 0.005
viewer.settings['scale.reactions'] = 0.005/3
viewer.settings['scale.reactions'] = 0.005/6
# viewer.settings['scale.reactions'] = 0.005*2
viewer.settings['opacity.shapes'] =  0.3

viewer.draw_shape()
viewer.show()

