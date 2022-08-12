from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer

file_formdiagram = '/Users/mricardo/compas_dev/me/anagni/revision/form_min_good_sol.json'

form = FormDiagram.from_json(file_formdiagram)

t_w = form.thrust()/form.lumped_swt()
print('SWT:', form.lumped_swt())
print('T/W', t_w)

viewer = Viewer(form)
viewer.draw_thrust()
viewer.draw_shape()
viewer.draw_reactions()
viewer.draw_cracks()
viewer.show()

# form.to_json('/Users/mricardo/compas_dev/me/anagni/revision/form_min_good_sol.json')
