from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form

file_address = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form = FormDiagram.from_json(file_address)

for key in form.vertices():
    pz = form.vertex_attribute(key, 'pz')
    print(pz)

for key in form.vertices_where({'is_fixed': True}):
    print(key)
    pz = form.vertex_attribute(key, 'pz')
    print(pz)

plot_form(form, show_q=False, heights=True).show()
