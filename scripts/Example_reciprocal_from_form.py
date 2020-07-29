from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.viewers.thrust import view_thrust
import compas_tno

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax

import os

type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'
objective = 'max'

file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ objective + '.json')
form = FormDiagram.from_json(file_address)

# plot_form(form).show()
print(form)
form.overview_forces()

form, force = form.reciprocal_from_form(plot=False)
form.overview_forces()
# plot_form(form).show()
# form.plot()
# view_thrust(form).show()

file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ objective + '_force.json')
force.to_json(file_address)

# file_address = '/Users/mricardo/compas_dev/me/reformulation/orthogonal.json'
# form.to_json(file_address)
