from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters.form import plot_form_semicirculararch_xz
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
import os

type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [20, 16]
thickness_scaled_for = 'f'

folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# for thk, obj in [[20.454691527171835, 'min'], [49.99999999999999, 'min'], [49.99999999999999, 'max']]:
#     address = os.path.join(folder, title) + '_' + obj + '_thk_' + str(thk) + '.json'
#     print(address)
#     thk = thk/100
#     type_structure = 'dome'
#     discretisation = [20, 16]
#     n = 1
#     radius = 5.0

#     data_shape = {
#         'type': type_structure,
#         'thk': thk,
#         'discretisation': [discretisation[0]*n, discretisation[1]*n],
#         'center': [5.0, 5.0],
#         'radius': radius,
#         't': 0.0,
#         'expanded': True
#     }

#     dome = Shape.from_library(data_shape)
#     form = FormDiagram.from_json(address)
#     form.envelope_from_shape(dome)
#     form.overview_forces()
#     # form.to_json(address)

#     # address_plot = os.path.join(folder, title) + '_' + obj + '_thk_' + str(thk) + '_plot_' + thickness_scaled_for + '.pdf'
#     # plot_form(form, show_q=False, thick=thickness_scaled_for, cracks=True, save=address_plot).show()
#     # plot_form(form, show_q=False, thick='q', cracks=True).show()

#     # Plot the cross section

#     # invert sign of loads and forces

#     for u, v in form.edges():
#         q = form.edge_attribute((u, v), 'q')
#         form.edge_attribute((u, v), 'q', -q)

#     for key in form.vertices():
#         pz = form.vertex_attribute(key, 'pz')
#         form.vertex_attribute(key, 'pz', -pz)

#     compute_reactions(form)

#     form.attributes['Re'] = radius + thk/2
#     form.attributes['Ri'] = radius - thk/2
#     tol = 10e-3
#     address_plot_section = os.path.join(folder, title) + '_' + obj + '_thk_' + str(thk) + '_plot_' + 'section' + '.pdf'
#     plot_form_semicirculararch_xz(form, radius=0.06, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True, yrange=[radius-tol, radius+tol], save=address_plot_section).show()

#     from compas_tno.viewers import draw_thrust_as_lines
#     draw_thrust_as_lines(form).show()

address = '/Users/mricardo/compas_dev/me/general_opt/dome/radial_fd/mov_c_0.1/dome_radial_fd_discr_[20, 16]_t_thk_10.77604794596367.json'
form = FormDiagram.from_json(address)
thk = form.attributes['thk']
type_structure = 'dome'
discretisation = [20, 16]
n = 1
radius = 5.0
form.overview_forces()

compute_reactions(form)

from compas_tno.viewers import draw_thrust_as_lines
draw_thrust_as_lines(form).show()

tol = 10e-3
plot_form_semicirculararch_xz(form, radius=0.06, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True, yrange=[radius-tol, radius+tol]).show()
