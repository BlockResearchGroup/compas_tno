
import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_normals
from compas_tno.viewers import view_mesh
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.datastructures import MeshDos
from compas.datastructures import mesh_delete_duplicate_vertices
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from numpy import array
from scipy import rand

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

sols = {}
# for discretisation in [14]:
# for A in [1.00, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.20]:
for A in [1.05]:

    # Basic parameters

    discretisation = 14
    span_x = 10.0
    span_y = 10.0
    span = max(span_x, span_y)
    xy_span = [[0, span_x], [0, span_y]]
    xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
    thk = 0.5
    error = 0.0
    k = 1.0
    n = 1
    type_structure = 'crossvault'
    type_formdiagram = 'cross_fd'
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
    # zt = crossvault_middle_update(xyz[:, 0], xyz[:, 1],  t,  xy_span=xy_span_shape)
    zub, zlb = crossvault_ub_lb_update(xyz[:, 0], xyz[:, 1], thk, t,  xy_span=xy_span_shape)
    intrados = MeshDos.from_mesh(form)
    i = 0
    for key in intrados.vertices():
        intrados.vertex_attribute(key, 'z', zlb[i])
        i += 1

    intrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
    intrados.store_normals(correct_creases=True)

    extrados = intrados.offset_mesh(n=thk, direction='up')
    middle = intrados.offset_mesh(n=thk/2, direction='up')

    vault = Shape.from_meshes(intrados, extrados, middle=middle, data={'type': 'general', 't': 0.0, 'thk': thk})

    area = vault.middle.area()
    swt = vault.compute_selfweight()

    print('Interpolated Volume Data:')
    print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
    print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

    ze = extrados.vertices_attributes('z')
    print('min ze:', min(ze))

    # view_shapes(vault).show()
    # view_mesh(vault.intrados, normals=True).show()

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
    optimiser.data['variables'] = ['ind', 'zb', 't']
    optimiser.data['objective'] = 't'
    optimiser.data['thickness_type'] = 'intrados'  # 'variable', 'intrados'
    optimiser.data['min_thk'] = 0.0
    optimiser.data['max_thk'] = thk*1.5
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
        thk_min = optimiser.fopt
        data_shape['thk'] = thk_min
        # vault = Shape.from_library(data_shape)
        # form.envelope_from_shape(vault)

        # sols[str(discretisation)] = thk_min
        sols[str(A)] = thk_min

        # print('Solved:', discretisation, thk_min)
        print('Solved:', A, thk_min)

        # plot_form(form, show_q=False, simple=True, cracks=True).show()
        # view_solution(form, vault).show()

    else:
        # print('Not Solved:', discretisation)
        print('Not Solved:', A)

print(sols)

for sol in sols:
    # discr = int(sol)
    discr = float(sol)
    print('{0}, {1}'.format(discr, sols[sol]))
