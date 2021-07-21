
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
from compas_tno.viewers import view_meshes
from compas_tno.shapes import MeshDos
from compas.datastructures import mesh_delete_duplicate_vertices
from compas.geometry import normalize_vector
from compas.geometry import sum_vectors
from compas.geometry import norm_vector
from compas.geometry import scale_vector
from compas.geometry import angle_vectors
from compas_tno.plotters import save_pointcloud
from scipy import rand
import math

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

sols = {}
# for discretisation in [10, 12, 14, 16, 18, 20]:
for discretisation in [10]:

    # Basic parameters

    thk = 0.5
    error = 0.0
    span = 10.0
    xc = yc = span/2
    k = 1.0
    ro = 1.0
    n = 2
    type_structure = 'crossvault'
    type_formdiagram = 'cross_fd'
    gradients = True  # False

    # ----------------------- Shape Analytical ---------------------------

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation*n,
        'xy_span': [[0, span], [0, k*span]],
        't': 0.0,
        'expanded': False,
    }

    analytical_shape = Shape.from_library(data_shape)
    analytical_shape.ro = ro

    area_analytical = analytical_shape.middle.area()
    swt_analytical = analytical_shape.compute_selfweight()

    print('Analytical Self-weight is:', swt_analytical)
    print('Analytical Area is:', area_analytical)

    # ----------------------- Point Cloud -----------------------

    xy = []
    points_ub = []
    points_lb = []

    for i in range(n * discretisation + 1):
        for j in range(n * discretisation + 1):
            xy.append([i * span / (n * discretisation), j * span / (n * discretisation)])

    z_ub = analytical_shape.get_ub_pattern(xy).reshape(-1, 1)
    z_lb = analytical_shape.get_lb_pattern(xy).reshape(-1, 1)

    for i in range(len(xy)):
        points_lb.append([xy[i][0], xy[i][1], float(z_lb[i])])
        points_ub.append([xy[i][0], xy[i][1], float(z_ub[i])])

    # title = 'square_crossvault'
    # json = '/Users/mricardo/compas_dev/me/min_thk/pointcloud/' + title + '.json'
    # save_pointcloud(points_lb, points_ub, json)

    # triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
    # view_shapes_pointcloud(triangulated_shape).show()

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
    vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub, middle=analytical_shape.middle, data={'type': 'general', 't': 0.0, 'thk': thk})
    vault.store_normals()
    vault.ro = ro

    vault.intrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
    vault.extrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])

    vault.intrados.store_normals(correct_creases=True)
    vault.extrados.store_normals(correct_creases=True)

    for key in vault.intrados.vertices_where({'is_crease': True}):
        x, y, _ = vault.intrados.vertex_coordinates(key)
        n = vault.intrados.vertex_attribute(key, 'n')
        print('INTRA: x, y, n, norm', x, y, n, norm_vector(n))
    for key in vault.extrados.vertices_where({'is_crease': True}):
        x, y, _ = vault.extrados.vertex_coordinates(key)
        n = vault.extrados.vertex_attribute(key, 'n')
        print('EXTRA x, y, n, norm', x, y, n, norm_vector(n))

    from compas_plotters import MeshPlotter
    plotter = MeshPlotter(vault.intrados)
    plotter.draw_edges()
    plotter.draw_vertices(text={key: norm_vector(vault.intrados.vertex_attribute(key, 'n')) for key in vault.intrados.vertices()})
    plotter.show()

    plotter = MeshPlotter(vault.extrados)
    plotter.draw_edges()
    plotter.draw_vertices(text={key: norm_vector(vault.extrados.vertex_attribute(key, 'n')) for key in vault.extrados.vertices()})
    plotter.show()

    # view_normals(vault).show()

    area = vault.middle.area()
    swt = vault.compute_selfweight()

    print('Interpolated Volume Data:')
    print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
    print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

    form.selfweight_from_shape(vault)
    # form.selfweight_from_shape(analytical_shape)

    # --------------------- 3. Create Starting point with TNA ---------------------

    # form = form.initialise_tna(plot=False)
    form.initialise_loadpath()
    # plot_form(form).show()

    # --------------------- 4. Create Minimisation Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.data['library'] = 'Scipy'
    optimiser.data['solver'] = 'SLSQP'
    optimiser.data['constraints'] = ['funicular', 'envelope']
    optimiser.data['variables'] = ['ind', 'zb', 'n']
    optimiser.data['objective'] = 'n'
    optimiser.data['printout'] = True
    optimiser.data['plot'] = False
    optimiser.data['find_inds'] = True
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

    if optimiser.exitflag == 0:
        n_reduction = -1 * analysis.optimiser.fopt
        thk_min = thk - 2 * n_reduction
        data_shape['thk'] = thk_min
        vault = Shape.from_library(data_shape)
        form.envelope_from_shape(vault)

        folder = os.path.join('/Users/mricardo/compas_dev/me', 'max_n', type_structure, type_formdiagram)
        os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_offset-method'
        save_form = os.path.join(folder, title)

        form.to_json(save_form + '_min_thk_' + optimiser.data['objective'] + '_' + str(thk_min) + '.json')

        sols[str(discretisation)] = thk_min

        print('Solved:', discretisation, thk_min)

    else:

        print('Not Solved:', discretisation)

print(sols)

for sol in sols:
    discr = int(sol)
    print('{0}, {1}'.format(discr, sols[sol]))
