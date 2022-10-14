import os
os.environ["USE_PROPACK"] = "1"

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import equilibrium_fdm
from compas_plotters import Plotter
from compas.datastructures import Mesh
from compas_tno.utilities import apply_bounds_on_q
from compas.colors import Color
from compas.geometry import Translation
from compas.geometry import Scale
from compas.geometry import Point
from compas.geometry import distance_point_point_xy
from compas_view2.shapes import Arrow
import compas_tno
from numpy import zeros
import os

from compas.datastructures import mesh_weld

# Geometry parameters

radius = 5.0
thk = 0.5
discretisation = [8, 10]
discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

xp, yp = 7.5, 5.0

# Parameters Optimisations

obj = 'max_load'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']  # , 'displ_map'
variables = ['q', 'zb', 'lambdv']  # zb, lambdv
features = ['fixed']
starting_point = 'loadpath'
make_video = False
autodiff = False

# Create shape/diagram

dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.5)
folder = '/Users/mricardo/compas_dev/compas_tno/data/'
folder = '/Users/mricardo/compas_dev/me/pattern/singular/dome/'
folder = '/Users/mricardo/compas_dev/me/max_load/dome/apex/dome/'

for prob in ['A2-sym']:  # ['A1', 'B1', 'C1', 'D1', 'E1']  # Then H and then > D3-diag whithout independent edges (overnight) >>> H2 BEST
    # Try this script and maybe try also with independent edges...

    mesh_file = folder + 'mesh-' + prob + '.json'

    mesh = Mesh.from_json(mesh_file)
    print('mesh faces:', mesh.number_of_faces())

    mesh_weld(mesh)

    form = FormDiagram.from_mesh(mesh)
    print('form faces:', form.number_of_faces())

    form.delete_boundary_edges()
    form.set_boundary_supports()

    # form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)  # , diagonal=True, partial_diagonal='straight')

    # # helper to select node to load

    # form = FormDiagram.create_circular_spiral_form(discretisation=discretisation, radius=radius)

    # Maximum load magnitude

    max_load_mult = 600.0
    n = form.number_of_vertices()
    pzv = zeros((n, 1))

    bbox = form.bounding_box_xy()
    xmin, ymin = min([m[0] for m in bbox]), min([m[1] for m in bbox])

    if abs(xmin) > 10e-3 or abs(ymin) > 10e-3:
        dx, dy = -xmin, -ymin
        translation = Translation.from_vector([dx, dy, 0.0])
        form.transform(translation)

        bbox = form.bounding_box_xy()
        xmax, ymax = max([m[0] for m in bbox]), max([m[1] for m in bbox])
        print(bbox)
        print(xmax, ymax)

        scale = Scale.from_factors([10.0/xmax, 10.0/ymax, 0.0])
        form.transform(scale)

        bbox = form.bounding_box_xy()
        print(bbox)

    dist = 0.01
    for key in form.vertices():
        coords = form.vertex_coordinates(key)
        print(key, coords, distance_point_point_xy(coords, [xp, yp]))
        if distance_point_point_xy(coords, [xp, yp]) < dist:
            load_node = key

    if prob == 'D2-mod':
        load_node = 101

    # plotter = Plotter()
    # plotter.fontsize = 9
    # artist = plotter.add(form)
    # artist.draw_vertexlabels()
    # plotter.zoom_extents()
    # plotter.show()

    pzv[load_node] = -1.0

    # plotter = TNOPlotter(form)
    # plotter.draw_form(scale_width=False)
    # plotter.draw_supports()
    # plotter.add(Point(*form.vertex_coordinates(load_node)), facecolor=Color.red())
    # plotter.show()

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
    optimiser.settings['find_inds'] = False
    optimiser.settings['plot'] = False
    optimiser.settings['save_iterations'] = make_video
    optimiser.settings['solver_convex'] = 'CVXPY'
    optimiser.settings['max_iter'] = 1000
    optimiser.settings['autodiff'] = autodiff

    optimiser.settings['max_lambd'] = max_load_mult
    optimiser.settings['load_direction'] = pzv

    # Create analysis

    apply_bounds_on_q(form, qmin=-1000)

    analysis = Analysis.from_elements(dome, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    # analysis.apply_envelope_on_xy(c=-.1)
    analysis.apply_reaction_bounds()

    pz0 = form.vertex_attribute(0, 'pz')

    pzt = 0
    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        pzt += pz

    print('Total load of:', pzt)

    analysis.set_up_optimiser()

    ad = '/Users/mricardo/compas_dev/compas_tno/data/analysis-dome' + prob + '-with_inds.json'
    analysis.to_json(ad)

    # to view starting point
    # view = Viewer(form, dome)
    # view.draw_thrust()
    # view.draw_shape()
    # view.draw_force()
    # view.show()

    # to view starting point
    # view = Viewer(form, dome)
    # view.settings['camera.show.grid'] = False
    # view.settings['camera.distance'] = 35
    # view.settings['camera.target'] = [5, 5, 0]
    # view.settings['camera.rz'] = 0
    # view.draw_thrust()
    # view.draw_shape()
    # view.show()

    # plotter = TNOPlotter(form, figsize=(14, 6))
    # plotter.settings['size.edge.max_thickness'] = 8.0
    # plotter.draw_form(scale_width=True)
    # plotter.draw_supports()
    # plotter.draw_force()
    # plotter.show()

    analysis.run()

    fopt = analysis.optimiser.fopt
    exitflag = analysis.optimiser.exitflag

    print('Exitflag is:', exitflag)

    pc = fopt/pzt

    print('Percentage of load added is:', round(pc*100, 3), '%')

    lastanalysis = compas_tno.get('form.json')
    form.to_json(lastanalysis)

    if exitflag == 0:

        foldersave = os.path.join(folder, 'mesh-' + prob)
        os.makedirs(foldersave, exist_ok=True)
        title = 'dome' + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '_pct_stw_' + str(pc) + '.json'
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

    plotter = TNOPlotter(form, shape=dome)
    plotter.draw_form()
    plotter.draw_supports()
    plotter.draw_cracks()
    plotter.show()

    # to view end point
    length = 2.0
    x, y, z = form.vertex_coordinates(load_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    view = Viewer(form, dome)

    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['camera.rz'] = -180
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
        animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, shape=dome, interval=150, record=save_gif, jump_each=interv)
