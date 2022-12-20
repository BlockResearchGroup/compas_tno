from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.problems import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_shape
from compas.colors import Color
import math

# form = FormDiagram.create_circular_radial_form(discretisation=[16, 20], diagonal=True)
# form = FormDiagram.create_circular_spiral_form(discretisation=[16, 20])

# form = FormDiagram.create_fan_form(discretisation=discr)
# form = FormDiagram.create_cross_form(discretisation=discr)
# form = FormDiagram.create_cross_with_diagonal(discretisation=14, fix='all')
# form = FormDiagram.create_cross_diagonal(discretisation=14, fix='all')
# form = FormDiagram.create_parametric_form(discretisation=14, lambd=0.3, fix='all')
# form = FormDiagram.create_cross_form(discretisation=14, fix='all')
# form = FormDiagram.create_delta_form(discretisation=14, delta=0.5)

for lambd in [0, 1]:

    # form = FormDiagram.create_parametric_form(discretisation=14, lambd=lambd)

    form = FormDiagram.create_delta_form(discretisation=14, delta=0.714)

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_cracks()
    plotter.draw_supports(color=Color.red(), size=8)
    plotter.show()
