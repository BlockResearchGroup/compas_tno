from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.problems import initialise_form

form = FormDiagram.create_cross_form(discretisation=14, fix='all')

initialise_form(form)

inds = list(form.edges_where({'is_ind': True}))
print(len(inds))

plt = TNOPlotter(form)
plt.draw_form_independents()
plt.draw_supports()
plt.show()
