from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.viewers import Viewer
from compas_tno.utilities import update_json
from compas_tno.plotters import TNOPlotter
from compas_plotters import Plotter

ad = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk.json'

ad = update_json(ad)

form = FormDiagram.from_json(ad)
force = ForceDiagram.from_formdiagram(form)

force = reciprocal_from_form(form)

ad2 = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_force.json'
force.to_json(ad2)

# plotter = Plotter()
# plotter.add(form)
# plotter.show()

# plotter = Plotter()
# plotter.add(force)
# plotter.show()

# form.plot()
# force.plot()

# print('Continue code 1')

# plotter = TNOPlotter(form)
# plotter.show_solution()

# print('Continue code 2')

# plotter = MeshPlotter(form, figsize=(8, 8))
# plotter.draw_edges(keys=[edge for edge in form.edges_where({'_is_edge': True})], width=1.0)
# plotter.draw_vertices(keys=form.fixed(), facecolor='FF0000')
# plotter.show()

# plotter = MeshPlotter(force, figsize=(8, 8))
# plotter.draw_edges(keys=[edge for edge in force.edges()], width=1.0)
# plotter.show()

# plot_form(form, show_q=True).show()

# for u, v in form.edges():
#     qi = form.edge_attribute((u, v), 'q')
#     fi = form.edge_attribute((u, v), 'q')
#     print(fi)
#     # print(qi)

# viewer = Viewer(form)
# viewer.settings['size.edge.max_thickness'] = viewer.settings['size.edge.max_thickness']/4.0
# viewer.draw_thrust()
# viewer.view_cracks()
# viewer.show()
