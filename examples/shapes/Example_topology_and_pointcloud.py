from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.viewers import view_shapes

from scipy import rand

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
span = 10.0
k = 1.0
n = 2
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
gradients = True

# ----------------------- Shape Analytical ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

analytical_shape = Shape.from_library(data_shape)

area_analytical = analytical_shape.middle.area()
swt_analytical = analytical_shape.compute_selfweight()

print('Analytical Self-weight is:', swt_analytical)
print('Analytical Area is:', area_analytical)

# ----------------------- Populate Point Cloud ---------------------------

xy = []
points_ub = []
points_lb = []
error = 0.10

for i in range(n * discretisation + 1):
    for j in range(n * discretisation + 1):
        xy.append([i * span / (n * discretisation), j * span / (n * discretisation)])

from compas_tno.utilities import get_shape_ub_pattern
from compas_tno.utilities import get_shape_lb_pattern
z_ub = get_shape_ub_pattern(analytical_shape, xy).reshape(-1, 1) + error * (2 * rand(len(xy), 1) - 1)
z_lb = get_shape_lb_pattern(analytical_shape, xy).reshape(-1, 1) + error * (2 * rand(len(xy), 1) - 1)

for i in range(len(xy)):
    points_lb.append([xy[i][0], xy[i][1], float(z_lb[i])])
    points_ub.append([xy[i][0], xy[i][1], float(z_ub[i])])

# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a point cloud

vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub)

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

view_shapes(vault).show()
