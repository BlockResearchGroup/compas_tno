
import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_normals
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.datastructures import MeshDos
from compas.datastructures import mesh_delete_duplicate_vertices
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from numpy import array
from scipy import rand
import math

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

disc = []
discretisations = [12]
degs = [30]

for discretisation in discretisations:
    sols = {}

    type_formdiagram = 'fan_fd'
    type_structure = 'crossvault'
    span = span_x = span_y = 10.0
    xy_span = [[0, span_x], [0, span_y]]
    k = 1.0
    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners',
    }

    form = FormDiagram.from_library(data_diagram)
    print('Form Diagram Created!')
    j = 0

    for deg in degs:

        A = 1/math.cos(math.radians(deg))
        print('A, deg, discretisation:', A, deg, discretisation)

        xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
        thk = 0.40
        n = 1
        type_structure = 'crossvault'
        gradients = True  # False
        t = 0.0

        # ----------------------- Shape Analytical ---------------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': discretisation*n,
            'xy_span': xy_span_shape,
            't': t,
        }

        analytical_shape = Shape.from_library(data_shape)
        area_analytical = analytical_shape.middle.area()
        swt_analytical = analytical_shape.compute_selfweight()
        print('Analytical Self-weight is:', swt_analytical)
        print('Analytical Area is:', area_analytical)

        # ----------------------- Form Diagram ---------------------------

        data_diagram = {
            'type': type_formdiagram,
            'xy_span': [[0, span], [0, k*span]],
            'discretisation': discretisation,
            'fix': 'corners',
        }

        form = FormDiagram.from_library(data_diagram)
        print('Form Diagram Created!')
        # plot_form(form, show_q=False, fix_width=False).show()

        # ------- Create shape given a topology and a point cloud --------

        # roots - not considering real middle
        # vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub)
        # improved, considers the real middle

        xyz = array(form.vertices_attributes('xyz'))
        zt = crossvault_middle_update(xyz[:, 0], xyz[:, 1],  t,  xy_span=xy_span_shape)
        zub, zlb = crossvault_ub_lb_update(xyz[:, 0], xyz[:, 1], thk, t,  xy_span=xy_span_shape)
        intrados = MeshDos.from_mesh(form)
        extrados = MeshDos.from_mesh(form)
        middle = MeshDos.from_mesh(form)
        i = 0
        for key in middle.vertices():
            intrados.vertex_attribute(key, 'z', zlb[i])
            extrados.vertex_attribute(key, 'z', zub[i])
            middle.vertex_attribute(key, 'z', zt[i])
            i += 1

        vault = Shape.from_meshes(intrados, extrados, middle=middle, data={'type': 'general', 't': 0.0, 'thk': thk})

        vault.intrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
        vault.extrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])

        vault.intrados.store_normals(correct_creases=True)
        vault.extrados.store_normals(correct_creases=True)

        area = vault.middle.area()
        swt = vault.compute_selfweight()

        print('Interpolated Volume Data:')
        print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
        print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

        # view_shapes(vault).show()

        form.selfweight_from_shape(vault)
        # form.selfweight_from_shape(analytical_shape)

        # SAVE FOLDER

        folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', 'numerical_n', type_structure, type_formdiagram)
        os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_deg=' + str(deg)
        save_form = os.path.join(folder, title)

        # --------------------- 3. Create Starting point with TNA ---------------------

        # form = form.initialise_tna(plot=False)
        # if j == 0:
        try:
            form = FormDiagram.from_json(save_form + '_lp_' + '.json')
        except:
            form.initialise_loadpath()
            form.to_json(save_form + '_lp_' + '.json')
        # plot_form(form).show()

        # --------------------- 4. Create Minimisation Optimiser ---------------------

        optimiser = Optimiser()
        optimiser.data['library'] = 'Scipy'
        optimiser.data['solver'] = 'SLSQP'
        optimiser.data['constraints'] = ['funicular', 'envelope']
        optimiser.data['variables'] = ['ind', 'zb', 'n']
        optimiser.data['objective'] = 'n'
        optimiser.data['thickness_type'] = 'constant'
        optimiser.data['min_thk'] = -thk
        optimiser.data['max_thk'] = 0.0
        optimiser.data['printout'] = True
        optimiser.data['plot'] = False
        optimiser.data['find_inds'] = True
        optimiser.data['max_iter'] = 5000
        optimiser.data['qmax'] = 1000.0
        optimiser.data['gradient'] = gradients
        optimiser.data['jacobian'] = gradients

        # --------------------- 5. Set up and run analysis ---------------------

        analysis = Analysis.from_elements(vault, form, optimiser)
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()
        analysis.run()

        j += 1

        if optimiser.exitflag == 0:
            n_reduction = -1 * analysis.optimiser.fopt
            thk_min = thk - 2 * n_reduction
            data_shape['thk'] = thk_min

            form.to_json(save_form + '_min_thk_' + optimiser.data['objective'] + '_' + str(thk_min) + '.json')

            # sols[str(discretisation)] = thk_min
            # sols[str(A)] = thk_min
            sols[str(deg)] = thk_min

            # print('Solved:', discretisation, thk_min)
            # print('Solved:', A, thk_min)
            print('Solved:', deg, thk_min)

            # plot_form(form, show_q=False, simple=True, cracks=True).show()
            # view_solution(form, vault).show()

        else:
            # print('Not Solved:', discretisation)
            # print('Not Solved:', A)
            print('Not Solved:', deg)

    print(type_formdiagram)
    print(discretisation)
    print(sols)
    for sol in sols:
        discr = float(sol)
        print('{0}, {1}'.format(discr, sols[sol]))
    disc.append(sols)

i = 0
for sols in disc:
    print('discretisation:', discretisations[i])
    for sol in sols:
        discr = float(sol)
        print('{0}, {1}'.format(discr, sols[sol]))
    i += 1