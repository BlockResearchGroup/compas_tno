from compas_tno.diagrams import FormDiagram
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas_tno.plotters import plot_form
from compas_tno.viewers import view_thrust
from compas.datastructures import mesh_delete_duplicate_vertices

mesh = Mesh.from_json('/Users/mricardo/compas_dev/compas_tno/data/dome/singular/test.json')
mesh_delete_duplicate_vertices(mesh)
form = FormDiagram.from_mesh(mesh)
print(form)
form.set_boundary_supports()
form.delete_boundary_edges()

plot_form(form, show_q=False, fix_width=True).show()

form = form.initialise_tna(plot=True, alpha=90.0)
plot_form(form, show_q=False).show()

form.to_json('/Users/mricardo/compas_dev/compas_tno/data/dome/singular/form1.json')
view_thrust(form).show()
