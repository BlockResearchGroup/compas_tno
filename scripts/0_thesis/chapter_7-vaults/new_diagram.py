from compas_tno.diagrams import FormDiagram
from compas_tno.problems import initialise_form
from compas_tno.utilities.form import slide_pattern_inwards
from compas_tno.algorithms import apply_sag
from compas_tno.plotters import TNOPlotter

form = FormDiagram.create_cross_form(discretisation=14)

apply_sag(form)

plot = TNOPlotter(form)
plot.draw_form(scale_width=False)
plot.show()

# initialise_form(form)

# path = '/Users/mricardo/compas_dev/me/pattern/new_parametric/form_{}.json'.format()

