from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.problems.initialize import initialize_tna
from compas_tno.problems.initialize import initialize_fdm
from compas_tno.utilities import apply_selfweight_from_thrust
from compas_tno.problems import initialize_loadpath
from compas.colors import Color
from compas_tno.viewers import Viewer

form : FormDiagram = FormDiagram.create_ortho_form(discretisation=6, fix='all')
apply_selfweight_from_thrust(form, thickness=1.0, density=20.0)

print('Form weight is:', form.lumped_swt())

plt = TNOPlotter(form)
plt.draw_form(scale_width=False, color=Color.black())
plt.draw_supports(color=Color.red())
plt.show()

# initialize_loadpath(form)
initialize_tna(form)

for edge in form.edges_where({'_is_edge': True}):
    q = form.edge_attribute(edge, 'q')
    form.edge_attribute(edge, 'q', 2*q)

initialize_fdm(form)

view : Viewer = Viewer(form, show_grid=False)
view.draw_thrust()
view.draw_reactions()
view.show()

path = '/Users/mricardo/compas_dev/me/thesis/networks/initial_network.json'
form.to_json(path)
