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
from compas_tno.viewers import view_mesh
from compas_tno.shapes import MeshDos
from scipy import interpolate
import json
from scipy import rand
from numpy import array

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.60
error = 0.0
span = 10.0
k = 1.0
n = 2
ro = 1.0
type_formdiagram = 'fan_fd'
discretisation = 10
gradients = True  # False

# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    # 'fix': 'all',
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
# plot_form(form, show_q=False, fix_width=False).show()

# ----------------------- Point Cloud -----------------------

# files = ['corners1', 'corners2', 'corners3']
files = ['corners1']
# files = ['continuous1', 'continuous2', 'continuous3']
# files = ['continuous3']

for file_name in files:
    folder = '/Users/mricardo/compas_dev/me/freeform/IASS/'
    pointcloud = folder + file_name + '.json'

    tol = 10e-4

    middle_pts = []
    xy_middle = []
    z_middle = []

    with open(pointcloud) as json_file:
        data = json.load(json_file)
        for key, pt in data['target'].items():
            xy_middle.append(pt[:2])
            z_middle.append(pt[2])
            middle_pts.append(pt)

    # triangulated_shape = Shape.from_middle_pointcloud(middle_pts, thk=thk)
    # view_normals(triangulated_shape).show()
    # view_shapes_pointcloud(triangulated_shape).show()

    # ------- Create shape given a topology and a point cloud --------

    form_xyz = array(form.vertices_attributes('xyz'))
    zt = interpolate.griddata(xy_middle, z_middle, form_xyz[:, :2], method='linear')
    middle = MeshDos.from_mesh(form)
    i = 0
    for key in middle.vertices():
        middle.vertex_attribute(key, 'z', zt[i])
        i += 1

    vault = Shape.from_middle(middle, thk=thk, treat_creases=False)
    vault.ro = ro

    # view_shapes_pointcloud(vault).show()
    # view_normals(vault).show()

    area = vault.middle.area()
    swt = vault.compute_selfweight()

    print('Interpolated Volume Data:')
    print('Self-weight is: {0:.2f}'.format(swt))
    print('Area is: {0:.2f}'.format(area))

    # view_shapes(vault).show()

    form.selfweight_from_shape(vault)
    # form.selfweight_from_shape(analytical_shape)

    # --------------------- 3. Create Starting point with TNA ---------------------

    # form = form.form_update_with_parallelisation(plot=False)
    form.initialise_loadpath()
    # plot_form(form).show()

    # --------------------- 4. Create Minimisation Optimiser ---------------------

    optimiser = Optimiser()
    optimiser.settings['library'] = 'Scipy'
    optimiser.settings['solver'] = 'SLSQP'
    optimiser.settings['constraints'] = ['funicular', 'envelope']
    optimiser.settings['variables'] = ['ind', 'zb', 't']
    optimiser.settings['objective'] = 't'
    optimiser.settings['printout'] = True
    optimiser.settings['plot'] = False
    optimiser.settings['find_inds'] = True
    optimiser.settings['qmax'] = 1000.0
    optimiser.settings['gradient'] = gradients
    optimiser.settings['jacobian'] = gradients

    # --------------------- 5. Set up and run analysis ---------------------

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    analysis.run()

    for key in form.vertices_where({'is_fixed': True}):
        z = form.vertex_attribute(key, 'z')
        ub = form.vertex_attribute(key, 'ub')
        lb = form.vertex_attribute(key, 'lb')
        print('key z, ub, lb:', key, z, ub, lb)

    if optimiser.exitflag == 0:

        thk_min = form.attributes['thk']
        plot_form(form, show_q=False, simple=True, cracks=True).show()

        folder_save = os.path.join(folder, type_formdiagram)
        os.makedirs(folder_save, exist_ok=True)
        title = file_name + '_' + type_formdiagram + '_discr_' + str(discretisation)
        save_form = os.path.join(folder_save, title)
        # form.to_json(save_form + '_min_thk_' + optimiser.settings['objective'] + '_' + str(thk_min) + '.json')

        view_solution(form, vault).show()
