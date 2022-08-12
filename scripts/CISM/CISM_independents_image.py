from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
from compas_tno.problems import initialise_form
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.utilities import apply_envelope_from_shape
from compas.colors import Color

delta = 1.0
span = 10.0
xspan = yspan = [0.0, span]
xspan_vault = yspan_vault = [- delta, span + delta]

form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM.json')
# form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM-2.json')
# form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM-3.json')

compute_reactions(form)

vault = Shape.create_crossvault(xy_span=[xspan_vault, yspan_vault])
apply_envelope_from_shape(form, vault)

vault = None
vault = Shape.from_formdiagram_and_attributes(form)

# M = initialise_form(form)
# adapt_problem_to_fixed_diagram(M, form, printout=True)

# plotter = TNOPlotter(form)
# plotter.settings['size.vertex'] = 8.0
# plotter.settings['color.edges.independent'] = Color.blue()
# plotter.draw_form_independents()
# plotter.draw_supports(color=Color.red())
# plotter.show()

# plotter = TNOPlotter(form)
# plotter.settings['size.vertex'] = 8.0
# plotter.settings['color.edges.independent'] = Color.blue()
# plotter.draw_force(show_independents=True)
# plotter.show()

viewer = Viewer(form, shape=vault, show_grid=False)
viewer.settings['camera.target'] = [5, 5, 0]
viewer.settings['camera.distance'] = 35
viewer.settings['scale.reactions'] = 0.005
viewer.settings['scale.reactions'] = 0.005/2
viewer.settings['scale.reactions'] = 0.005*2
viewer.settings['opacity.shapes'] =  0.3
viewer.draw_thrust()
viewer.draw_shape()
viewer.draw_reactions()
viewer.show()
