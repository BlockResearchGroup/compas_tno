from compas.datastructures import Mesh
from compas.datastructures import mesh_bounding_box_xy
from compas_plotters import MeshPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import project_mesh_to_middle


span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'fan_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 2 * discretisation

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 0.0,
}

vault = Shape.from_library(data_shape)

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)

# dual = form.build_dual()
dual = form.build_complete_dual()
project_mesh_to_middle(dual, vault)

# plotter = MeshPlotter(dual)
# plotter.draw_edges()
# plotter.draw_faces()
# plotter.draw_vertices(text={key: key for key in dual.vertices()})
# plotter.show()

view = Viewer(form, vault)
view.view_thrust()
# view.view_mesh()
view.view_mesh(dual)
# view.view_shape()
view.show()
