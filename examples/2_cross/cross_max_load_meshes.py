from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import equilibrium_fdm
from compas_tno.algorithms import apply_sag
from compas_plotters import Plotter
from compas.datastructures import Mesh
from compas_tno.utilities import apply_bounds_on_q
from compas.colors import Color
from compas.geometry import Translation
from compas.geometry import Scale
from compas.geometry import distance_point_point_xy
from compas_view2.shapes import Arrow
from numpy import zeros
import os

# Geometry parameters

delta = 1.0
span = 10.0
xspan = yspan = [0.0, span]
corners = [[0.0, 0.0], [span, 0.0], [span, span], [0.0, span]]
xspan_vault = yspan_vault = [- delta, span + delta]
thk = 0.5
discretisation = 20

# Parameters Optimisations

obj = 'max_load'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 'lambdv']  # zb, lambdv
features = ['fixed']
axis_sym = None  # [[[0, 5.0], [10, 5.0]], [[5.0, 0], [5.0, 10.0]], [[0.0, 0.0], [10.0, 10.0]]]
starting_point = 'loadpath'
make_video = False
autodiff = False
multiplier = -1.0

# Create shape/diagram

vault = Shape.create_crossvault(xy_span=[xspan_vault, yspan_vault], discretisation=discretisation, thk=thk)
folder = '/Users/mricardo/compas_dev/compas_tno/data/'
folder = '/Users/mricardo/compas_dev/me/pattern/singular/dome/'
folder = '/Users/mricardo/compas_dev/me/pattern/singular/crossvault/'

for prob in ['E2-web']:  # ['A1', 'B1', 'C1', 'D1', 'E1']

    mesh_file = folder + 'mesh-' + prob + '.json'

    mesh = Mesh.from_json(mesh_file)
    print('mesh faces:', mesh.number_of_faces())

    form = FormDiagram.from_mesh(mesh)
    print('form faces:', form.number_of_faces())
    print('form edges:', form.number_of_edges())

    # form.delete_boundary_edges()
    # form.set_boundary_supports()

    # form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)  # , diagonal=True, partial_diagonal='straight')

    # helper to select node to load
    plotter = Plotter()
    plotter.fontsize = 12
    artist = plotter.add(form)
    artist.draw_vertexlabels()
    plotter.zoom_extents()
    plotter.show()

    # form = FormDiagram.create_circular_spiral_form(discretisation=discretisation, radius=radius)

    # Maximum load magnitude

    max_load_mult = 1500.0
    n = form.number_of_vertices()
    pzv = zeros((n, 1))

    if prob == 'A1':  # robin patterns
        load_node = 66
    elif prob == 'A2':  # robin patterns
        load_node = 138
    elif prob == 'A2-web':  # robin patterns
        load_node = 272
    elif prob == 'B2':  # robin patterns
        load_node = 12
    elif prob == 'B2-web':  # robin patterns
        load_node = 188
    elif prob == 'C2':  # robin patterns
        load_node = 172
    elif prob == 'D2-web':  # robin patterns
        load_node = 229
    elif prob == 'E2-web':  # robin patterns
        load_node = 229

    pzv[load_node] = multiplier

    move_pattern_to_origin(corners, form)

    apply_sag(form, boundary_force=20.0)

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_supports()
    plotter.show()

    # Create optimiser

    optimiser = Optimiser()
    optimiser.settings['objective'] = obj
    optimiser.settings['solver'] = solver
    optimiser.settings['constraints'] = constraints
    optimiser.settings['variables'] = variables
    optimiser.settings['features'] = features
    optimiser.settings['starting_point'] = starting_point
    optimiser.settings['derivative_test'] = False
    optimiser.settings['printout'] = True
    optimiser.settings['plot'] = True
    optimiser.settings['save_iterations'] = make_video
    optimiser.settings['max_iter'] = 1500
    optimiser.settings['axis_sym'] = axis_sym
    optimiser.settings['autodiff'] = autodiff

    optimiser.settings['max_lambd'] = max_load_mult
    optimiser.settings['load_direction'] = pzv

    # Create analysis

    apply_bounds_on_q(form, qmin=-1000)

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    # analysis.apply_envelope_on_xy(c=-.1)

    pz0 = form.vertex_attribute(0, 'pz')

    pzt = 0
    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        pzt += pz

    plotter = Plotter()
    plotter.fontsize = 9
    artist = plotter.add(form)
    artist.draw_vertexlabels(text={key: round(form.vertex_attribute(key, 'pz'), 1) for key in form.vertices()})
    plotter.zoom_extents()
    plotter.show()

    print('Total load of:', pzt)

    vault0 = Shape.from_formdiagram_and_attributes(form)

    analysis.set_up_optimiser()

    # 2D plot form diagram
    plotter = TNOPlotter(form, figsize=(14, 6))
    plotter.settings['size.edge.max_thickness'] = 8.0
    plotter.draw_form(scale_width=True)
    plotter.draw_supports()
    plotter.draw_force()
    plotter.show()

    analysis.run()

    fopt = analysis.optimiser.fopt
    exitflag = analysis.optimiser.exitflag

    print('Exitflag is:', exitflag)

    pc = abs(multiplier * fopt / pzt)

    print('Selfweight is:', abs(pzt))
    print('Percentage of load added is:', round(pc*100, 3), '%')

    if exitflag == 0:

        foldersave = os.path.join(folder, 'mesh-' + prob)
        os.makedirs(foldersave, exist_ok=True)
        title = 'vault' + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '_pct_stw_' + str(pc) + '.json'
        save_form = os.path.join(foldersave, title)
        form.to_json(save_form)

        print('Solution Saved at:', save_form)

    # view = Viewer(form, dome)
    # view.draw_thrust()
    # view.draw_shape()
    # view.draw_force()
    # view.draw_cracks()
    # view.draw_reactions()
    # view.show()

    plotter = TNOPlotter(form, figsize=(14, 6))
    plotter.settings['size.edge.max_thickness'] = 8.0
    plotter.draw_form(scale_width=True)
    plotter.draw_supports()
    plotter.draw_cracks()
    plotter.draw_force()
    plotter.show()

    plotter = TNOPlotter(form, shape=vault)
    plotter.draw_form()
    plotter.draw_supports()
    plotter.draw_cracks()
    plotter.show()

    # to view end point
    length = 2.0
    x, y, z = form.vertex_coordinates(load_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    view = Viewer(form, vault0)

    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['scale.reactions'] = 0.001
    view.draw_thrust()
    view.draw_shape()
    view.draw_cracks()
    view.app.add(arrow, color=(0, 0, 0))
    view.draw_reactions()
    view.show()

    if make_video:

        from compas_tno.viewers import animation_from_optimisation
        from compas_tno.algorithms import reciprocal_from_form
        import compas_tno

        DATA_XFORM = compas_tno.get('Xform.json')
        DATA_XFORCE = compas_tno.get('Xforce.json')

        force = reciprocal_from_form(form)

        save_gif = folder + '-iterations.gif'
        interv = 1
        if len(DATA_XFORM) > 100:
            interv = round(len(DATA_XFORM)/50)
        animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, shape=vault, interval=150, record=save_gif, jump_each=interv)
