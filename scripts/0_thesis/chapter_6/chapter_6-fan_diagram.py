from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.problems import initialise_form

form = FormDiagram.create_fan_form(discretisation=14)

initialise_form(form)

plt = TNOPlotter(form)
plt.draw_form_independents()
plt.draw_supports()
plt.show()
