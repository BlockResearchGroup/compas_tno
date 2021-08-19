
from compas_tno.diagrams import FormDiagram
from compas.datastructures import Mesh
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_thrust_as_lines
from compas_tno.plotters import plot_form
import compas
import compas_rv2

formjson = '/Users/mricardo/compas_dev/me/freeform/tom/form2.json'
# targetjson = '/Users/mricardo/compas_dev/me/freeform/tom/target.json'

form = FormDiagram.from_json(formjson)
# target = Mesh.from_json(targetjson)

# vertices, faces = form.to_vertices_and_faces()
# form = FormDiagram.from_vertices_and_faces(vertices, faces)

n = form.number_of_vertices()

print('-'*20)
print('Number of edges (form):', form.number_of_edges())
print('Number of vertices (form):', form.number_of_vertices())

pzt = 0.0
thk = 1.0
# for key in target.vertices():
#     zt = target.vertex_attribute(key, 'z')
#     form.vertex_attribute(key, 'z', zt)
#     area = form.vertex_area(key)
#     pz = thk * area
#     pzt += pz
#     form.vertex_attribute(key, 'pz', -pz)

count_fixed = 0
thk = 2.0
for key in form.vertices():
    zt = form.vertex_attribute(key, 'z')
    form.vertex_attribute(key, 'target', zt)
    area = form.vertex_area(key)
    pz = thk * area
    pzt += pz
    form.vertex_attribute(key, 'pz', -pz)
    if form.vertex_attribute(key, 'is_anchor') is True:
        form.vertex_attribute(key, 'is_fixed', True)
        count_fixed += 1

print('Fixed Nodes:', count_fixed)
print('Free Nodes:', n - count_fixed)

for u, v in form.edges():
    q = form.edge_attribute((u, v), 'q')
    form.edge_attribute((u, v), 'q', -q)

print('Total load applied:', pzt)
print('-'*20)

# formtno = '/Users/mricardo/compas_dev/me/freeform/tom/form_tno.json'
# form.to_json(formtno)

plot_form(form, show_q=False, max_width=1.0).show()

view_thrust_as_lines(form).show()
# view_thrust(target).show()
view_thrust(form).show()
