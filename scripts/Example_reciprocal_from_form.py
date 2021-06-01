from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.viewers.thrust import view_thrust
from compas_tno.viewers import view_shapes
import compas_tno

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax

import os

type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'
objective = 'min'

# file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ objective + '.json')
file_address = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=6.4843/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_15.207373123795708.json'

file_address = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=6.1147/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_17.723600813371625.json'

file_address = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.24_discr_[8, 20]_min.json'
# file_address = '/Users/mricardo/compas_dev/me/shape_comparison/dome/radial_spaced_fd/dome_radial_spaced_fd_discr_[8, 20]_min_thk_50.0.json'
# file_address = '/Users/mricardo/compas_dev/me/shape_comparison/dome/radial_spaced_fd/dome_radial_spaced_fd_discr_[8, 20]_max_thk_45.0.json'

form = FormDiagram.from_json(file_address)

plot_form(form).show()
plot_form(form, show_q=False, cracks=True).show()
print(form)
form.overview_forces()


corners = list(form.vertices_where({'is_fixed': True}))
form.vertices_attribute('is_anchor', True, keys=corners)
react = {'_rx': {}, '_ry': {}}
for key in form.vertices_where({'is_fixed': True}):
    react['_rx'][key] = form.vertex_attribute(key, '_rx')
    react['_ry'][key] = form.vertex_attribute(key, '_ry')

print('num. edges:', form.number_of_edges())
print('num. faces:', form.number_of_faces())
print('num. vertices:', form.number_of_vertices())
form.update_boundaries()
print('num. edges:', form.number_of_edges())
print('num. faces:', form.number_of_faces())
print('num. vertices:', form.number_of_vertices())

# from compas_plotters import MeshPlotter
# plotter = MeshPlotter(form)
# # plotter.draw_edges(text={edge: round(form.edge_length(*edge)*form.edge_attribute(edge, 'q'), 1) for edge in form.edges()})
# plotter.draw_edges()
# plotter.draw_vertices()
# # plotter.draw_faces()
# plotter.draw_faces(text={key: key for key in form.faces()})
# plotter.show()

from compas_tno.algorithms import force_update_from_form

from compas_tno.diagrams import ForceDiagram
force = ForceDiagram.from_formdiagram(form)

# print('Plot of Dual')
# force.plot()

# force_update_from_form(force, form)

form, force = form.reciprocal_from_form(plot=False)
form.overview_forces()

print('Plot of Dual')
force.plot()

# from compas_plotters import MeshPlotter
# plotter = MeshPlotter(force)
# plotter.draw_edges(text={edge: round(force.edge_length(*edge), 1) for edge in force.edges()})
# plotter.draw_vertices()
# # plotter.draw_faces()
# plotter.show()

max_length = 0
for edge in force.edges():
    length = force.edge_length(*edge)
    if length > max_length:
        max_length = length

print('max_length', max_length)   # Fix this such rhar

# form, force = form.reciprocal_from_form(plot=False)
# form.overview_forces()

# print('Plot of Dual - TNA')
# force.plot()

# max_length = 0
# for edge in force.edges():
#     length = force.edge_length(*edge)
#     if length > max_length:
#         max_length = length

# print('max_length', max_length)

faces_unloaded = []

for face in form.faces_where({'is_loaded': False}):
    faces_unloaded.append(face)
print('faces unloaded:', faces_unloaded)
for face in faces_unloaded:
    form.delete_face(face)

print('num. edges:', form.number_of_edges())
print('num. faces:', form.number_of_faces())
print('num. vertices:', form.number_of_vertices())

dxy = 1.0
from compas_plotters import MeshPlotter
plotter = MeshPlotter(force)
# plotter.draw_edges(text={edge: round(force.edge_length(*edge), 1) for edge in force.edges()})
plotter.draw_edges()
plotter.draw_vertices(text={key: key for key in force.vertices()})
# plotter.draw_faces()
plotter.show()

# Find out the nodes that are superimposed, and move them inwards... Then update the form diagram

# plot_form(form).show()
# form.plot()
view_thrust(form).show()
# view_shapes(form).show()

# file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ objective + '_force.json')
# force.to_json(file_address)

# file_address = '/Users/mricardo/compas_dev/me/reformulation/orthogonal.json'
# form.to_json(file_address)
