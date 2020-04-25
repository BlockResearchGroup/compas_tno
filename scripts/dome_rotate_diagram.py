from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_shapes
import compas
from compas.geometry import matrix_from_axis_and_angle
from compas.datastructures import Mesh
from compas.datastructures import mesh_transform
from math import pi

# T = matrix_from_axis_and_angle([0, 1, 0], pi / 4)

from compas.geometry import Frame
from compas.geometry import Transformation
from compas.geometry import Sphere
from compas.geometry import Point
sphere = Sphere(Point(1, 1, 1), 5)
lambda_mult = 0.20

frame = Frame([5, 5, 0], [1, 0, -lambda_mult], [0, 1, 0])
T = Transformation.from_frame(frame)
sphere_transformed = sphere.transformed(T)
print(sphere)
# view_thrust(sphere).show()

data_diagram = {
    'type': 'radial_spaced_fd',
    'D': 3.0,
    'center': [5.0, 5.0],
    'radius': 5.0,
    'discretisation': [10, 20],
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)

data_shape = {
    'type': 'dome_polar',
    'thk': 0.15,
    'discretisation': [10, 20],
    'center': [5.0, 5.0],
    'radius': 5.0,
    't' : 0.0
}

shape = Shape.from_library(data_shape)
# view_shapes(shape).show()

analysis = Analysis.from_form_and_shape(form, shape)
analysis.apply_target()
analysis.apply_envelope()

for key in form.vertices():
    form.vertex_attribute(key, 'z', form.vertex_attribute(key, 'target'))

# view_thrust(form).show()

from compas.datastructures import mesh_smooth_area
from compas.datastructures import mesh_smooth_centroid
fixed = [key for key in form.vertices() if form.vertex_attribute(key, 'is_fixed') == True]
# mesh_smooth_area(form, fixed=fixed)
mesh_smooth_centroid(form, fixed=fixed)

view_thrust(form).show()

tmesh = form.copy()
# mesh_transform(tmesh, T)
plot_form(tmesh, show_q=False).show()
view_thrust(tmesh).show()
