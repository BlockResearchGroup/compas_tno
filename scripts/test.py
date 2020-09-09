
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
address = '/Users/mricardo/compas_dev/me/shape_comparison/domicalvault/cross_fd/domicalvault_cross_fd_discr_10_min_thk_50.0.json'

form = FormDiagram.from_json(address)
plot_form(form, show_q=True).show()

print('PZ')
for key in form.vertices():
    print(form.vertex_attribute(key, 'pz'))

print('Q')
for u, v in form.edges():
    print(form.edge_attribute((u,v), ' q'))
