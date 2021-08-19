
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form

# -----------------------------------------------------------
# --------- CREATING CROSS W DIAGONAL FORM DIAGRAM ------------
# -----------------------------------------------------------

data = {
    'type': 'cross_diagonal',
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
print(form)
plot_form(form, show_q=False).show()
