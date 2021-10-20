
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_normals
from compas_tno.viewers import view_mesh
from compas_tno.viewers import Viewer
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import Viewer
from compas_tno.shapes import MeshDos
from compas.geometry import normalize_vector
from compas.geometry import sum_vectors
from compas.geometry import norm_vector
from compas.geometry import scale_vector
from compas.geometry import angle_vectors
from compas.geometry import centroid_points
from compas_tno.shapes.crossvault import crossvault_middle_update
import math
from numpy import array
from copy import deepcopy
import os

# ----------------------------------------------------------------------
# --------- NUMERICAL MIN THK WITH SPRING ANGLE - CROSSVAULT -----------
# ----------------------------------------------------------------------

# A = 1.1  # Parameter similar to a spring angle

disc = []
# discretisations = [10, 12, 14, 16, 18, 20]
discretisations = [16]
degs = [0, 5, 10, 15, 20, 25, 30, 35, 40]

for discretisation in discretisations:
    sols = {}

    type_formdiagram = 'cross_fd'
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

    # for A in [1.00, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.20]:
    # for A in [1.00, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.20]:
    for deg in degs:
    # for A in [1.00]:

        A = 1/math.cos(math.radians(deg))
        print('A, deg, discretisation:', A, deg, discretisation)
        thk = 0.5

        xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
        k = 1.0
        n = 2
        type_structure = 'crossvault'
        type_formdiagram = 'cross_fd'
        gradients = True
        t = 0.0

        # --------------------
        # Shape
        # --------------------

        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': discretisation*n,
            'xy_span': xy_span_shape,
            't': t,
        }

        analytical_shape = Shape.from_library(data_shape)
        swt = analytical_shape.compute_selfweight()
        print('Selfweight computed:', swt)
        print('Vault geometry created!')

        # --------------------
        # PointCloud in the middle
        # --------------------

        # xy = []
        # points_central = []

        # for i in range(n * discretisation + 1):
        #     for j in range(n * discretisation + 1):
        #         xy.append([i * span / (n * discretisation), j * span / (n * discretisation)])

        # zt = analytical_shape.get_middle_pattern(xy).reshape(-1, 1)

        # for i in range(len(xy)):
        #     points_central.append([xy[i][0], xy[i][1], float(zt[i])])

        # vault = Shape.from_middle_pointcloud(points_central, topology=form, thk=thk, treat_creases=True)

        # --------------------
        # Mesh in the middle
        # --------------------

        xyz = array(form.vertices_attributes('xyz'))
        zt = crossvault_middle_update(xyz[:, 0], xyz[:, 1],  t,  xy_span=xy_span_shape)
        middle = MeshDos.from_mesh(form)
        i = 0
        for key in middle.vertices():
            middle.vertex_attribute(key, 'z', zt[i])
            i += 1

        vault = Shape.from_middle(middle, thk=thk, treat_creases=True)

        # --------------------
        # Mesh as a percentage
        # --------------------

        # from numpy import array
        # XY = array(vault.middle.vertices_attributes('xyz'))[:, :2]
        # zub = analytical_shape.get_ub_pattern(XY)
        # zlb = analytical_shape.get_lb_pattern(XY)

        # vault.middle.scale_normals_with_ub_lb(zub, zlb)
        # vault.extrados, vault.intrados = vault.middle.offset_up_and_down(n=1.0)

        # vault.middle.scale_normals
        # view_shapes(vault).show()
        # view_shapes(analytical_shape).show()

        # --------------------
        # Initialise
        # --------------------

        # form = form.form_update_with_parallelisation(plot=False)
        form.selfweight_from_shape(vault)
        form.envelope_from_shape(vault)
        if j == 0:
            form.initialise_loadpath()
        # plot_form(form).show()

        # --------------------
        # Optimiser
        # --------------------

        optimiser = Optimiser()
        optimiser.settings['library'] = 'Scipy'
        optimiser.settings['solver'] = 'slsqp'
        optimiser.settings['constraints'] = ['funicular', 'envelope']
        # optimiser.settings['constraints'] = ['funicular', 'envelope', 'symmetry']
        optimiser.settings['variables'] = ['ind', 'zb', 't']
        optimiser.settings['objective'] = 't'
        optimiser.settings['thickness_type'] = 'constant'
        optimiser.settings['min_thk'] = 0.0
        optimiser.settings['max_thk'] = thk*2.0
        optimiser.settings['printout'] = True
        optimiser.settings['plot'] = False
        optimiser.settings['find_inds'] = True
        optimiser.settings['qmax'] = 1000.0
        optimiser.settings['gradient'] = gradients
        optimiser.settings['jacobian'] = gradients

        # --------------------
        # Analysis
        # --------------------

        analysis = Analysis.from_elements(vault, form, optimiser)
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.set_up_optimiser()
        analysis.run()

        j+=1

        if optimiser.exitflag == 0:

            thk_min = form.attributes['thk']

            folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', 'numerical_t', type_structure, type_formdiagram)
            os.makedirs(folder, exist_ok=True)
            title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_deg=' + str(deg)
            save_form = os.path.join(folder, title)

            form.to_json(save_form + '_min_thk_' + optimiser.settings['objective'] + '_' + str(thk_min) + '.json')

            # sols[str(discretisation)] = thk_min
            # sols[str(A)] = thk_min
            sols[str(deg)] = thk_min

            # print('Solved:', discretisation, thk_min)
            # print('Solved:', A, thk_min)
            print('Solved:', deg, thk_min)

            # plot_form(form, show_q=False, simple=True, cracks=True).show()
            # view_solution(form, vault).show()
            # view_shapes(vault).show()
            # view_normals(vault).show()
        else:
            # print('Not Solved:', discretisation)
            # print('Not Solved:', A)
            print('Not Solved:', deg)

    print(discretisation)
    print(sols)
    disc.append(sols)

print(disc)

i = 0
for sols in disc:
    print('discretisation:', discretisations[i])
    for sol in sols:
        discr = float(sol)
        print('{0}, {1}'.format(discr, sols[sol]))
    i += 1

