from compas_tno.algorithms.graphstatics import form_update_with_parallelisation
from compas_tno.diagrams import FormDiagram
from compas_tno.problems import initialise_form
from compas_tno.algorithms import form_update_with_parallelisation
from compas_tno.plotters import TNOPlotter


discr = 4
x0 = y0 = 0.0
x1 = y1 = 10.0
delta = 10.0/discr/2

# XXXX

form = FormDiagram.create_cross_form(discretisation=discr, fix='corners')

M = initialise_form(form, printout=True)
print('Shape of Eq. Matrix:', M.E.shape)

force = form_update_with_parallelisation(form, printout=True)

plotter = TNOPlotter(form, force=force, figsize=(16, 6))
plotter.draw_form_independents()
plotter.draw_force()
plotter.show()




# # XXXX

# form = FormDiagram.create_cross_form(discretisation=discr, fix='corners')

# slide_diagram(delta, y0, y1, form)

# M = initialise_form(form, printout=True)
# print('Shape of Eq. Matrix:', M.E.shape)

# force = form_update_with_parallelisation(form, printout=True)

# plotter = TNOPlotter(form, force=force, figsize=(16, 6))
# plotter.draw_form_independents()
# plotter.draw_force()
# plotter.show()

# # XXXX

# form = FormDiagram.create_cross_with_diagonal(discretisation=discr, fix='corners')

# slide_diagram(delta, y0, y1, form)

# M = initialise_form(form, printout=True)
# print('Shape of Eq. Matrix:', M.E.shape)

# force = form_update_with_parallelisation(form, printout=True)

# plotter = TNOPlotter(form, force=force, figsize=(16, 6))
# plotter.draw_form_independents()
# plotter.draw_force()
# plotter.show()



