from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import view_mesh

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
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
