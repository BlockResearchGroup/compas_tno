from compas_tno.diagrams import FormDiagram
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.algorithms import apply_sag
from compas_tno.plotters import TNOPlotter
from compas.colors import Color

from compas_tno.utilities.form import slide_diagram
import math

data = {
    'type': 'cross_fd',  # 'cross_fd', 'ortho_fd'
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 12,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)

# plotter = TNOPlotter(form)
# plotter.draw_form(scale_width=False, color=Color.black())
# plotter.draw_supports(color=Color.red())
# plotter.show()

# shear = Shear.from_angle_direction_plane(angle=math.radians(30))
from compas.geometry import Shear
from compas.geometry import Scale
from compas.geometry import Frame
from compas.geometry import Point


# from compas.geometry import cross_vectors
# angle = 0.8
# direction = [0.1, 0.2, 0.3]
# point = [4, 3, 1]
# normal = cross_vectors(direction, [1, 0.3, -0.1])
# S = Shear.from_angle_direction_plane(angle, direction, (point, normal))
# form.transform(S)
# point = Point(2, 5, 0)
# frame = Frame(point, (1, 0, 0), (0, 1, 0))
# # points = [point, Point(2, 10, 0)]
# # S = Scale.from_factors([2.] * 3, frame)
# # frame = Frame(point, [1, 0, 0], [0, 1, 0])
# S = Scale.from_factors([1.5, 1.1, 0.2], frame)
# form.transform(S)

slide_diagram(form, delta=-0.5)

plotter = TNOPlotter(form)
# plotter.draw_mesh()
plotter.draw_form(scale_width=False, color=Color.black())
plotter.draw_supports(color=Color.red())
plotter.show()

for bf in [20.0]:

    form = FormDiagram.from_library(data)

    for key in form.vertices():
        form.vertex_attribute(key, 'pz', -1.0)

    apply_sag(form, boundary_force=bf)

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False, color=Color.black())
    plotter.draw_supports(color=Color.red())
    plotter.show()

    # force = reciprocal_from_form(form)
    # # force = ForceDiagram.from_formdiagram(form)

    # plotter = MeshPlotter(force, figsize=(8, 8))
    # plotter.draw_edges()
    # plotter.show()
