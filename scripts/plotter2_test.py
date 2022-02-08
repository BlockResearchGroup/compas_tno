from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities.symmetry import apply_radial_symmetry

# filepath = '/Users/mricardo/compas_dev/compas_tno/data/form-minthk.json'
filepath = '/Users/mricardo/compas_dev/compas_tno/data/form-compl.json'
form = FormDiagram.from_json(filepath)

apply_radial_symmetry(form)

filepath = '/Users/mricardo/compas_dev/compas_tno/data/shape.json'
shape = Shape.from_json(filepath)

# plotter = TNOPlotter(form, shape=shape)
plotter = TNOPlotter(form, figsize=(16, 8))

# plotter.draw_form_independents()
# plotter.draw_mesh()
plotter.draw_form_sym()
# plotter.draw_shape()
# plotter.draw_base_form()
# plotter.draw_form()
plotter.draw_supports()
# plotter.draw_cracks()
# plotter.draw_reactions()
# plotter.draw_force()
# plotter.highlight_vertices([50, 53], show_forcepolygon=True)
# plotter.draw_base_form()
# plotter.draw_form_xz()
# plotter.draw_shape_xz()
# plotter.zoom_extends()
plotter.show()
