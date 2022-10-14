from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.problems import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_thrust
from compas.colors import Color


# form_address = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_entropy.json'
form_address = '/Users/mricardo/compas_dev/me/freeform/tom/form_lpopt.json'

# form_address = '/Users/mricardo/compas_dev/me/freeform/IASS/continuous1.json'

form = FormDiagram.from_json(form_address)

form = FormDiagram.create_cross_form(discretisation=20, fix='all')

form = FormDiagram.create_fan_form(discretisation=20, fix='corners')

m = form.number_of_real_edges()
n = form.number_of_supports()

print('Edges, supports', m, n)

# form.boundary
# form.delete_boundary_edges()

# plt = TNOPlotter(form)
# plt.draw_form(scale_width=False, color=Color.black())
# plt.draw_supports(color=Color.red())
# plt.show()

# apply_selfweight_from_thrust(form, density=20.0)
# initialize_loadpath(form, printout=True, solver_convex='MATLAB')

# # save = '/Users/mricardo/compas_dev/me/freeform/tom/form_lpopt.json'
# # form.to_json(save)

# plt = TNOPlotter(form)
# plt.draw_form()
# plt.draw_supports(color=Color.red())
# plt.show()

# view = Viewer(form)
# view.settings['camera.target'] = [-7.5, 0.5, 0.0]
# view.draw_form()
# view.show()
