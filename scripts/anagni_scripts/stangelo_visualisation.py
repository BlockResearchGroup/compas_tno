from compas_view2.objects import Arrow
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
from compas_tno.problems import initialise_form
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.utilities import apply_envelope_from_shape
from compas.colors import Color
import compas_tno

path = compas_tno.open_dialog()

form = FormDiagram.from_json(path)
max_load = True

swt = form.lumped_swt()
thrust = form.thrust()
print('T/W:', round(thrust/swt, 2))

vault = Shape.from_formdiagram_and_attributes(form)

# # Plot form force
# plotter = TNOPlotter(form, figsize=(14, 6))
# plotter.settings['size.edge.max_thickness'] = 8.0
# plotter.settings['color.edges.form'] = Color.black()
# plotter.settings['color.vertex.supports'] = Color.red()
# plotter.draw_form(scale_width=False)
# plotter.draw_supports()
# plotter.draw_force()
# plotter.show()

# from compas_plotters import Plotter
# plotter = Plotter()
# artist = plotter.add(form)
# artist.draw_vertexlabels()
# plotter.show()

# # Classic plot
# plotter = TNOPlotter(form)
# plotter.show_solution()

viewer = Viewer(form, shape=vault, show_grid=False)
viewer.settings['camera.target'] = [2.8, 1.6, 12]
viewer.settings['camera.distance'] = 18
viewer.settings['scale.reactions'] = 0.005 * 4
viewer.settings['scale.loads'] = 0.005 * 4 * 10
viewer.settings['opacity.shapes'] = 0.3
viewer.settings['scale.edge.thk_absolute'] = 2.5
viewer.draw_thrust(absolute_scale=True)
viewer.draw_cracks()
viewer.draw_shape()
viewer.draw_reactions()

if max_load:
    length = 0.5
    xf = 5.744
    dist_load = 0.1
    xload = xf/3
    dict_loads = {}
    for i, key in enumerate(form.vertices()):
        xi, _, _ = form.vertex_coordinates(key)
        if abs(xi - xload) < dist_load:
            dict_loads[key] = key
    print(dict_loads)

    for key in dict_loads:
        x, y, z = form.vertex_coordinates(key)
        z += length + 0.2
        arrow = Arrow([x, y, z], [0, 0, -length])
        viewer.add(arrow, color=Color.black(), opacity=0.8)

# if objective == 'Ecomp-linear':
#     for key in dict_settlement:
#         x, y, z = form.vertex_coordinates(key)
#         dXbi = dict_settlement[key]
#         z -= 0.2
#         arrow = Arrow([x, y, z], dXbi)
#         if norm_vector(dXbi) > 0.01:
#             viewer.add(arrow, color=Color.black(), opacity=0.8)

viewer.show()


# from compas_view2.objects import Arrow
# loaded_node = 143
# length = 2.0
# x, y, z = form.vertex_coordinates(loaded_node)
# z += length + 0.1
# arrow = Arrow([x, y, z], [0, 0, -length])
# viewer.app.add(arrow, color=(0, 0, 0))

# viewer.draw_cracks()
# viewer.draw_shape()
# viewer.draw_reactions()
# viewer.show()

# Classic plot
plotter = TNOPlotter(form)
plotter.show_solution()
