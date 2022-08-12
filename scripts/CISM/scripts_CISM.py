from compas_tno.algorithms.graphstatics import form_update_with_parallelisation
from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.algorithms import form_update_with_parallelisation
from compas_tno.plotters import TNOPlotter
from compas.geometry import Vector
from compas.geometry import angle_vectors
from compas_plotters import Plotter
from compas.colors import Color
from compas.geometry import Translation
from compas_tno.viewers import Viewer

path = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM.json'
# path_force = '/Users/mricardo/compas_dev/compas_tno/data/CISM/force-CISM.json'

form = FormDiagram.from_json(path)
fixed = []
for key in form.vertices_where({'is_fixed': True}):
    form.vertex_attribute(key, 'z', 0.0)
    fixed.append(key)


ad = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/force-0.json'
force0 = FormDiagram.from_json(ad)

plotter = Plotter()
plotter.add(
    form,
    show_faces=False,
    show_vertices=True,
    vertices = fixed,
    edgecolor={edge: Color.grey() for edge in form.edges()},
    vertexcolor={key: Color.red() for key in fixed},
    edgewidth=1.0,
    vertexsize=10,
)
plotter.add(force0,
            show_faces=False,
            show_vertices=False,
            edgecolor={edge: Color.from_hex('#0092D2') for edge in force0.edges()},
            edgewidth=2.0
            )
plotter.zoom_extents()
plotter.show()

form.edges_attribute('is_ind', False)

angles = {}
for edge in form.edges():
    p1, p2 = form.edge_coordinates(*edge)
    ve = Vector(p2[0] - p1[0], p2[1] - p1[1], 0)
    v0 = Vector(1, 0, 0)
    angle = angle_vectors(ve, v0, deg=True)
    angles[edge] = round(angle)
    # print(angle)
    if abs(abs(angle) - 45) < 1 or abs(abs(angle) - 135) < 1:
        form.edge_attribute(edge, 'is_ind', True)

# plotter = TNOPlotter(form)
# plotter.draw_form_independents(show_text=False)
# plotter.formartist.draw_edgelabels(text=angles)
# plotter.show()

translation = [-10, 0, 0]
T = Translation.from_vector(translation)

# kmax = 100
# plotter = Plotter()
# for k in range(kmax + 1):
# # for k in [0, 1, 2, 3, 4, 5, 10, 20, 40, 60, 80, 100]:
#     # if not k % 2:
#     #     continue
#     ad = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/force-' + str(k) + '.json'
#     force = ForceDiagram.from_json(ad)
#     force.edges_attribute('_is_edge', True)
#     colors = {edge: Color.from_hex('#0092D2') for edge in force.edges()}
#     width = 0.1
#     if k == 0:
#         force.transform(T)
#     if k == 100:
#         width = 0.5
#         colors = {edge: Color.black() for edge in force.edges()}
#     plotter.add(force,
#                 show_faces=False,
#                 show_vertices=False,
#                 edgecolor=colors,
#                 edgewidth=width)

# plotter.zoom_extents()
# plotter.show()

force = form_update_with_parallelisation(form, kmax=100, callback=None)

view = Viewer(form, force=force, show_grid=False)
view.settings['size.edge.max_thickness'] = 30.0
view.draw_thrust()
# view.draw_force()
view.show()

ad = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/form-after-opt.json'
form.to_json(ad)

# print(force.vertex_coordinates(0))

# force.to_json(path_force)
# print('Force save:', path_force)
