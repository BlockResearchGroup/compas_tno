from compas_tno.diagrams import FormDiagram
from compas_plotters import Plotter, plotter
from compas.geometry import Line
# from matplotlib.collections import LineCollection
# from compas_plotters.artists import PlotterArtist
from compas_tno.plotters2 import FormPlotter
import compas_tno

json = compas_tno.get('form-test.json')

form = FormDiagram.from_json(json)

plotter = FormPlotter(form)
plotter.draw_form()
plotter.draw_cracks()
plotter.show()

# scale = 2.0

# plt = Plotter(figsize=(8, 8))

# max_q = abs(min(form.edges_attribute('q')))
# print(form.edges_attribute('q'))
# print(max_q)

# lines = []
# widths = []
# colors = []

# edgewidths = {}
# edgecolors = {}

# for edge in form.edges():
#     qi = abs(form.edge_attribute(edge, 'q'))
#     thk = qi / max_q * scale
#     p1, p2 = form.edge_coordinates(*edge)
#     line = Line(p1[:2], p2[:2])
#     # plt.add(line)
#     lines.append([p1[:2], p2[:2]])
#     widths.append(thk)
#     colors.append((1, 0, 0))
#     edgewidths[edge] = thk
#     edgecolors[edge] = (1, 0, 0)

# # collection = LineCollection(
# #             lines,
# #             linewidths=widths,
# #             colors=colors,
# #             linestyle='solid',
# #             alpha=1.0,
# #         )

# # plt.axes.add_collection(collection)

# meshartist = plt.add(form,
#                      vertices=form.fixed(),
#                      show_edges=True,
#                      show_faces=False,
#                      show_vertices=True,
#                      edgewidth=edgewidths,
#                      edge_color=edgecolors,
#                      vertexcolor=(0, 0, 0),
#                      vertexsize=10
#                      )

# # artist = PlotterArtist(form,
# #                      vertices=form.fixed(),
# #                      show_edges=True,
# #                      show_faces=False,
# #                      show_vertices=True,
# #                      edgewidth=edgewidths,
# #                      edge_color=edgecolors,
# #                      vertexcolor=(0, 0, 0),
# #                      vertexsize=10
# #                      )

# # plt._artists = [artist]
# # artist.draw_edges(color=edgecolors)

# plt.zoom_extents()
# plt.show()
