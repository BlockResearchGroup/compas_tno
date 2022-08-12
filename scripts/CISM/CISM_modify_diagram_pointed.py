from compas_tno.diagrams import FormDiagram
from compas_plotters import Plotter
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities import form_add_lines_support

form = FormDiagram.create_cross_form(discretisation=16)

yc = 5.0

plotter = Plotter()
artist = plotter.add(form)
artist.draw_vertexlabels()
plotter.show()

loaded_node = 142

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 154

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 162

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 170

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 170

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 176

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

loaded_node = 175

supports = [key for key in form.vertices_where({'is_fixed': True}) if form.vertex_coordinates(key)[1] < yc]

form, loaded_node = form_add_lines_support(form, loaded_node, supports)

plotter = Plotter()
artist = plotter.add(form)
artist.draw_vertexlabels()
plotter.show()


plotter = TNOPlotter(form)
plotter.draw_form()
plotter.draw_supports()
plotter.show()

form.to_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-directpath.json')
