
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.problems import initialise_form
from compas_tno.problems import initialise_problem
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

discretisation = 14
span = 10.0
type_structure = 'crossvault'
thk = 0.50
data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, span]],
    't': 0.0,
}

sols = {}
# for discretisation in [14]:

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
#     # plotter.save(save_img)
#     plotter.show()

#     # If required to find independent edges:

#     shape = Shape.from_library(data_shape)

#     form.selfweight_from_shape(shape)

#     initialise_form(form, printout=True)

#     title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_ind.pdf'
#     save_img = os.path.join(folder, title)
#     print(save_img)

#     plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
#     plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})], color={key: '0000FF' for key in form.edges_where({'is_ind': True})}, width={key: 2.5 for key in form.edges_where({'is_ind': True})})
#     plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.1, facecolor='FF0000')
#     # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': False})], radius=0.025, facecolor='000000')
#     plotter.save(save_img)
#     plotter.show()


### --------- amiens FORM DIAGRAM

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

    # plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
    # plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
    # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
    # plotter.save(save_img)
    # plotter.show()


    # If required to find independent edges:

    type_structure = 'crossvault'
    thk = 0.50
    data_shape['xy_span'] = [[0, span_y], [0, span_y]]

    shape = Shape.from_library(data_shape)
    form.selfweight_from_shape(shape)
    initialise_form(form, printout=True)

    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_ind.pdf'
    save_img = os.path.join(folder, title)
    print(save_img)

    plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
    plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})], color={key: '0000FF' for key in form.edges_where({'is_ind': True})}, width={key: 2.5 for key in form.edges_where({'is_ind': True})})
    plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='FF0000')
    # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': False})], radius=0.025, facecolor='000000')
    plotter.save(save_img)
    plotter.show()

