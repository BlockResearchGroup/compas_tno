from compas_tna.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import zlq_from_qid

file_address = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form = FormDiagram.from_json(file_address)

for key in form.vertices():
    form.vertex_attribute(key, 'pz', 10.0)
    pz = form.vertex_attribute(key, 'pz')
    print(pz)

plot_form(form).show()

from compas.geometry import mirror_points_line
