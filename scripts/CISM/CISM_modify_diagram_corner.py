from compas_tno.diagrams import FormDiagram
from compas_plotters import Plotter
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities import form_add_lines_support

# Quadrants map
# -----
# Q1 Q2
# Q3 Q4
# -----
#

# path = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-Ecomp-linear-16-thk-0.5-corner-diagonal.json'

form = FormDiagram.create_cross_form(discretisation=16)
# form = FormDiagram.from_json(path)

yc = 5.0
xc = 5.0

# plotter = Plotter()
# artist = plotter.add(form)
# artist.draw_vertexlabels()
# plotter.show()

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.show()

delta = 0.1

for key in form.vertices_where({'is_fixed': False}):
    x, y, _ = form.vertex_coordinates(key)
    if x > xc and y < yc:
        xn = x +
