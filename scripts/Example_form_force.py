from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.viewers.thrust import view_thrust

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax

data = {
    'type': 'cross_fd',
    'xy_span': [[0,10],[0,10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
# corners = list(form.vertices_where({'is_fixed': True}))
# form.vertices_attribute('is_anchor', True, keys=corners)
# form.edges_attribute('fmin', 0.0)
# form.update_boundaries(feet=2)
# form.plot()
# force = ForceDiagram.from_formdiagram(form)
# force.plot()
# horizontal_nodal(form, force, kmax=1000, display=False)
# force.plot()
# vertical_from_zmax(form, 5.0)
# # view_thrust(form).show()
# form = form.remove_feet()
# form.plot()
# view_thrust(form).show()

form = form.initialise_tna()
plot_form(form).show()
form.plot()
view_thrust(form).show()

file_address = '/Users/mricardo/compas_dev/me/reformulation/orthogonal.json'
form.to_json(file_address)
