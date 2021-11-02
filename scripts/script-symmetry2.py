from compas_tno.plotters import plot_symmetry
from compas_tno.plotters import plot_symmetry_vertices

from compas_tno.diagrams import FormDiagram

from compas_tno.utilities import apply_radial_symmetry
from compas_tno.utilities import apply_symmetry_from_axis

type_formdiagram = 'cross_fd'
span = 10.0
discretisation = 10

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)

hor_line = [[0.0, 5.0], [10.0, 5.0]]
ver_line = [[5.0, 0.0], [5.0, 10.0]]
diag_1 = [[0.0, 0.0], [10.0, 10.0]]
diag_2 = [[0.0, 10.0], [10.0, 0.0]]

# lines = [hor_line, ver_line]
lines = [hor_line, ver_line, diag_1]
# lines = [hor_line, hor_line]

# apply_symmetry(form, axis_symmetry=hor_line)
# apply_symmetry_from_axis(form, list_axis_symmetry=lines)
apply_radial_symmetry(form, list_axis_symmetry=lines)

plot_symmetry(form).show()
plot_symmetry_vertices(form).show()

