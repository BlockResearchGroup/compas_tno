from compas.utilities import i_to_red
from compas_plotters import MeshPlotter
import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_normals
from compas_tno.shapes import MeshDos
import math
from scipy import rand

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

sols = {}
for x_discr in [4, 8, 12, 16]:
    for y_discr in [12, 16, 20, 24, 28]:
        discretisation = [x_discr, y_discr]

        # Basic parameters

        thk = 0.5
        radius = 5.0
        type_structure = 'dome'
        type_formdiagram = 'radial_fd'
        # discretisation = [8, 20]
        center = [5.0, 5.0]
        gradients = False
        diagonal = False
        n = 1
        error = 0.00
        ro = 1.0

        # ----------------------- Shape Analytical ---------------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': [discretisation[0]*n, discretisation[1]*n],
            'center': center,
            'radius': radius,
            't': 0.0,
            'expanded': True
        }

        analytical_shape = Shape.from_library(data_shape)
        analytical_shape.ro = ro
        area_analytical = analytical_shape.middle.area()
        swt_analytical = analytical_shape.compute_selfweight()

        print('Analytical Self-weight is:', swt_analytical)
        print('Analytical Area is:', area_analytical)

        # analytical_shape.store_normals()
        # view_shapes(analytical_shape).show()

        # ----------------------- Form Diagram ---------------------------

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
        # print('Form Diagram Created!')
        # plot_form(form, show_q=False, fix_width=False).show()

        # ------- Create shape given a topology and a shape defined by meshes --------

        vault = Shape.from_meshes_and_formdiagram(form, analytical_shape.intrados, analytical_shape.extrados, middle=analytical_shape.middle,
                                                data={'type': 'general', 't': 0.0, 'thk': thk, 'base_structure': data_shape})
        vault.ro = ro

        area = vault.middle.area()
        swt = vault.compute_selfweight()

        print('Interpolated Volume Data:')
        print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
        print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

        vault.analytical_normals(assume_shape=data_shape)
        view_normals(vault).show()

        form.selfweight_from_shape(analytical_shape)
        form.envelope_from_shape(analytical_shape)
        # view_normals(vault).show()
        # view_shapes(vault).show()
        # view_shapes_pointcloud(vault).show()

        # --------------------- 3. Create Starting point with TNA ---------------------

        # form = form.initialise_tna(plot=False)
        form.initialise_loadpath()
        # plot_form(form).show()

        # --------------------- 4. Create Minimisation Optimiser ---------------------

        optimiser = Optimiser()
        optimiser.settings['library'] = 'Scipy'
        optimiser.settings['solver'] = 'SLSQP'
        optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
        optimiser.settings['variables'] = ['ind', 'zb', 'n']
        optimiser.settings['objective'] = 'n'
        optimiser.settings['printout'] = True
        optimiser.settings['plot'] = False
        optimiser.settings['find_inds'] = True
        optimiser.settings['qmax'] = 1000.0
        optimiser.settings['gradient'] = gradients
        optimiser.settings['jacobian'] = gradients

        # --------------------- 5. Set up and run analysis ---------------------

        analysis = Analysis.from_elements(vault, form, optimiser)
        # analysis.apply_selfweight()
        # analysis.apply_envelope()
        analysis.apply_reaction_bounds(assume_shape=data_shape)
        analysis.set_up_optimiser()
        analysis.run()

        n_reduction = - 1 * analysis.optimiser.fopt
        thk_min = thk - 2 * n_reduction
        print('Approx. Minimum THK:', thk_min)
        data_shape['thk'] = thk_min

        analytical_shape = Shape.from_library(data_shape)
        form.envelope_from_shape(analytical_shape)

        if optimiser.exitflag == 0:
            data_shape['thk'] = thk_min
            dome = Shape.from_library(data_shape)
            form.envelope_from_shape(dome)

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
            save_form = os.path.join(folder, title)

            form.to_json(save_form + '_min_thk_' + optimiser.settings['objective'] + '_' + str(thk_min) + '.json')

            sols[str(discretisation)] = thk_min

            print('Solved:', discretisation, thk_min)

        else:

            print('Not Solved:', discretisation)

        form.to_json(compas_tno.get('test.json'))


print(sols)
