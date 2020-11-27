
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from compas_tno.viewers import view_shapes
from compas_plotters import MeshPlotter
from copy import deepcopy
import os

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

# sols = {}
# for discretisation in [10, 12, 14, 16, 18, 20]:

#     thk = 0.5
#     span = 10.0
#     k = 1.0
#     n = 1
#     type_structure = 'crossvault'
#     type_formdiagram = 'fan_fd'
#     span_x = 10.0
#     span_y = 10.0

#     # ----------------------- Create Form Diagram ---------------------------

#     data_diagram = {
#         'type': type_formdiagram,
#         'xy_span': [[0, span_x], [0, span_y]],
#         'discretisation': discretisation,
#         'fix': 'corners',
#     }

#     form = FormDiagram.from_library(data_diagram)

#     # plot_form(form, simple=True, max_width=1, show_q=False).show()

#     folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'discretisations')
#     os.makedirs(folder, exist_ok=True)
#     title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '.pdf'
#     save_img = os.path.join(folder, title)

#     plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
#     plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
#     plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
#     plotter.save(save_img)
#     plotter.show()


sols = {}
for discretisation in [14]:

    thk = 0.5
    span = 10.0
    k = 1.0
    n = 1
    type_structure = 'crossvault'
    type_formdiagram = 'fan_fd'
    file_name = 'amiens_internet'
    span_x = 6.0
    span_y = 11.2

    # ----------------------- Create Form Diagram ---------------------------

    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span_x], [0, span_y]],
        'discretisation': discretisation,
        'fix': 'corners',
    }

    form = FormDiagram.from_library(data_diagram)

    # plot_form(form, simple=True, max_width=1, show_q=False).show()

    folder = os.path.join('/Users/mricardo/compas_dev/me', 'max_n', file_name, type_structure, type_formdiagram, 'min_max',  'discretisations')
    os.makedirs(folder, exist_ok=True)
    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_offset-method' + '.pdf'
    save_img = os.path.join(folder, title)

    plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
    plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
    plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
    plotter.save(save_img)
    plotter.show()

