from compas_tno.diagrams import FormDiagram
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas_tno.plotters import plot_form
from compas_tno.viewers import draw_thrust
from compas.datastructures import mesh_delete_duplicate_vertices


data = {'type': 'circular'}

form = FormDiagram.from_skeleton(data)
plot_form(form, show_q=False).show()
