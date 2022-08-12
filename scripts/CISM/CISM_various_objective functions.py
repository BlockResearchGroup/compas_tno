from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.algorithms import apply_sag
from compas_tno.utilities import form_add_lines_support
from compas.colors import Color
import compas_tno
from compas.geometry import Line
from compas_view2.shapes import Arrow
from compas.geometry import normalize_vector
from compas.geometry import norm_vector
from compas.geometry import Translation
from compas.geometry import Scale
from compas_plotters import Plotter
from numpy import array
from numpy import zeros
# from compas_tno.utilities import apply_envelope_from_shape

delta = 2.0
span = 10.0
xspan = yspan = [0.0, span]
xspan_vault = yspan_vault = [- delta, span + delta]
thk = 0.50
discr = 16
discr_shape = discr * 2
starting_point = 'loadpath'
sag = False
prob = None

# form = FormDiagram.create_cross_form(discretisation=discr)
# form = FormDiagram.create_cross_diagonal(discretisation=discr)
form = FormDiagram.create_cross_with_diagonal(discretisation=discr)
# form = FormDiagram.create_fan_form(discretisation=discr)


def displ_boundary(form, displ_type='corner-diagonal'):

    xspan, yspan = form.parameters['xy_span']
    Xc = [(xspan[1] - xspan[0])/2, (yspan[1] - yspan[0])/2, 0.0]
    print(Xc)

    vector_supports = []

    for key in form.vertices_where({'is_fixed': True}):
        x, y, _ = form.vertex_coordinates(key)

        dXbi = [0, 0, 0]

        if displ_type == 'outwards':
            sign = +1  # +1 for outwards / -1 for inwards
            dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), 0])  # 4 corners

        if displ_type == 'inwards':
            sign = -1  # +1 for outwards / -1 for inwards
            dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), 0])  # 4 corners

        if displ_type == 'corner-downwards':
            if x > Xc[0] and y < Xc[1]:             # vertical settlement
                dXbi = normalize_vector([0, 0, -1])

        if displ_type == 'side-downwards':
            if x > Xc[0]:             # vertical settlement
                dXbi = normalize_vector([0, 0, -1])

        if displ_type == 'corner-diagonal':
            if x > Xc[0] and y < Xc[1]:             # vertical settlement
                sign = +1
                dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), 0])

        vector_supports.append(dXbi)
        form.vertex_attribute(key, 'dXb', dXbi)

    dXb = array(vector_supports)
    print(dXb)

    return dXb


for thk in [0.50]:

    print('Optimisation for thk:', thk)

    vault = Shape.create_crossvault(xy_span=[xspan_vault, yspan_vault], discretisation=discr_shape, thk=thk)
    # form = FormDiagram.create_fan_form()

    objective = 'Ecomp-linear'  # try 'max' 'Ecomp-linear'
    solver = 'IPOPT'  # try SLSQP
    constraints = ['funicular', 'envelope']
    variables = ['q', 'zb']  # 'lambdv', 't'
    features = ['fixed']
    # axis_sym = [[[0.0, 5.0], [10.0, 5.0]], [[5.0, 0.0], [5.0, 10.0]], [[0.0, 0.0], [10.0, 10.0]]]
    optimiser = Optimiser()

    if objective == 'max_load':
        variables.append('lambdv')

        # form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-directpath.json')
        # loaded_node = 207  # pattern with multiple paths to support

        form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-appliedload.json')
        loaded_node = 241  # pattern with some paths to support

        # loaded_node = 143  # discretisation 16 center
        loaded_node = 139  # discretisation 16 web
        # loaded_node = 56  # discretisation 10 web

        # apply_sag(form, boundary_force=50)

        supports = [0, 271]

        plotter = Plotter()
        artist = plotter.add(form)
        artist.draw_vertexlabels()
        plotter.show()

        # form, loaded_node = form_add_lines_support(form, loaded_node, supports)

        # print('New loaded node:', loaded_node)

        # plotter = Plotter()
        # artist = plotter.add(form)
        # artist.draw_vertexlabels()
        # plotter.show()

        max_load_mult = 2000.0
        n = form.number_of_vertices()
        pzv = zeros((n, 1))
        pzv[loaded_node] = -1.0

        optimiser.settings['max_lambd'] = max_load_mult
        optimiser.settings['load_direction'] = pzv

    if objective == 'Ecomp-linear':

        form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM-3.json')
        # form = FormDiagram.from_json('//Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=6.4843/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_15.207373123795708.json')
        form.parameters['type'] = 'Braced'
        form.parameters['xy_span'] = [xspan, yspan]

        if sag:
            apply_sag(form, boundary_force=50)

        displ_type = 'corner-diagonal'
        dXb = displ_boundary(form, displ_type=displ_type)

        plotter = TNOPlotter(form)
        plotter.draw_form(scale_width=False)
        plotter.show()

        optimiser.settings['support_displacement'] = dXb

    optimiser.settings['library'] = solver
    optimiser.settings['solver'] = solver
    optimiser.settings['constraints'] = constraints
    optimiser.settings['variables'] = variables
    optimiser.settings['features'] = features
    optimiser.settings['objective'] = objective
    optimiser.settings['plot'] = False
    optimiser.settings['find_inds'] = False
    optimiser.settings['max_iter'] = 1500
    optimiser.settings['gradient'] = True
    optimiser.settings['jacobian'] = True
    optimiser.settings['printout'] = True
    optimiser.settings['starting_point'] = starting_point
    optimiser.settings['sym_loads'] = True
    # optimiser.settings['axis_sym'] = axis_sym
    optimiser.settings['max_thk'] = 0.1
    optimiser.settings['normalize_loads'] = False

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.set_up_optimiser()

    swt = form.lumped_swt()
    vault0 = Shape.from_formdiagram_and_attributes(form)

    analysis.run()

    if analysis.optimiser.exitflag == 0:
        path = compas_tno.get('')
        form_path = path + '/CISM/form-' + form.parameters['type'] + '-' + objective + '-' + str(discr) + '-thk-' + str(thk) + '.json'
        if sag:
            form_path = path + '/CISM/form-' + form.parameters['type'] + '-sag-' + objective + '-' + str(discr) + '-thk-' + str(thk)  + '.json'
        # form_path = path + '/CISM/form-temp.json'
        if objective == 'Ecomp-linear':
            form_path = path + '/CISM/form-' + form.parameters['type'] + '-' + objective + '-' + str(discr) + '-thk-' + str(thk) + '-' + displ_type + '.json'
            if sag:
                form_path = path + '/CISM/form-' + form.parameters['type'] + '-sag-' + objective + '-' + str(discr) + '-thk-' + str(thk) + '-' + displ_type + '.json'
        form.to_json(form_path)
        print('Form saved to:', form_path)
    else:
        optimiser.settings['starting_point'] = starting_point

    thrust = form.thrust()
    T_over_W = thrust/swt
    print('Thrust over Weight:', T_over_W)
    print('Form selfweight:', swt)

    if objective == 'max_load':
        fopt = optimiser.fopt
        P_over_W = abs(fopt/swt)
        print('Maximum Load over Weight:', P_over_W)

force = reciprocal_from_form(form)

# Plot form force
plotter = TNOPlotter(form, force=force, figsize=(14, 6))
plotter.settings['size.edge.max_thickness'] = 8.0
plotter.settings['color.edges.form'] = Color.black()
plotter.settings['color.vertex.supports'] = Color.red()
plotter.draw_form(scale_width=False)
plotter.draw_supports()
if force:
    plotter.draw_force()
plotter.show()

# Classic plot
plotter = TNOPlotter(form)
plotter.show_solution()

if objective == 't':
    vault0 = Shape.from_formdiagram_and_attributes(form)

viewer = Viewer(form, vault0)
viewer.settings['camera.show.grid'] = False
viewer.settings['camera.distance'] = 35
viewer.settings['camera.target'] = [5, 5, 0]
viewer.settings['scale.reactions'] = 0.005
viewer.settings['scale.reactions'] = 0.005/3
# viewer.settings['scale.reactions'] = 0.005*2
viewer.settings['opacity.shapes'] =  0.3
viewer.draw_thrust()

viewer.draw_thrust()
viewer.draw_cracks()
viewer.draw_shape()
viewer.draw_reactions()

if objective == 'max_load':
    length = 2.0
    x, y, z = form.vertex_coordinates(loaded_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    viewer.app.add(arrow, color=(0, 0, 0))

if objective == 'Ecomp-linear':
    for key in form.vertices_where({'is_fixed': True}):
        dXbi = form.vertex_attribute(key, 'dXb')
        x, y, _ = form.vertex_coordinates(key)
        z = form.vertex_attribute(key, 'lb') - 0.1
        if norm_vector(dXbi) > 0.01:
            arrow = Arrow([x, y, z], dXbi)
            viewer.app.add(arrow, color=(0, 0, 0))

viewer.show()
