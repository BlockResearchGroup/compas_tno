
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
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

sols = {}
for x_discr in [20]:  # More sensible  #[8, 12, 16, 20, 24] number parallel
    for y_discr in [16]:  # Less sensible  #[12, 16, 20, 24] number meridian
        discretisation = [x_discr, y_discr]

        thk = 0.5
        radius = 5.0
        type_structure = 'dome'
        type_formdiagram = 'radial_fd'
        # discretisation = [4, 12]
        ro = 1.0
        gradients = True
        n = 1

        # ----------------------- 1. Create Dome shape ---------------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': [discretisation[0]*n, discretisation[1]*n],
            'center': [5.0, 5.0],
            'radius': radius,
            't': 0.0,
            'expanded': True
        }

        dome = Shape.from_library(data_shape)
        dome.ro = ro
        swt = dome.compute_selfweight()
        print('Selfweight computed:', swt)
        print('Vault geometry created!')
        # view_shapes(dome).show()

        # ----------------------- 2. Create Form Diagram ---------------------------

        data_diagram = {
            'type': type_formdiagram,
            'center': [5.0, 5.0],
            'radius': radius,
            'discretisation': discretisation,
            'r_oculus': 0.0,
            'diagonal': False,
            'partial_diagonal': False,
        }

        form = FormDiagram.from_library(data_diagram)

        initialise_form(form, printout=True)

        # plot_form(form, simple=True, max_width=1, show_q=False).show()

        folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'discretisations')
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '.pdf'
        save_img = os.path.join(folder, title)

        plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
        plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
        plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='000000')
        # plotter.save(save_img)
        plotter.show()

        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_ind.pdf'
        save_img = os.path.join(folder, title)
        print(save_img)

        plotter = MeshPlotter(form, figsize=(10, 10), tight=True)
        plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})], color={key: '0000FF' for key in form.edges_where({'is_ind': True})}, width={key: 2.5 for key in form.edges_where({'is_ind': True})})
        plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.075, facecolor='FF0000')
        # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': False})], radius=0.025, facecolor='000000')
        plotter.save(save_img)
        plotter.show()
