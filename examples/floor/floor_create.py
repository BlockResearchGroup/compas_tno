from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import vertical_equilibrium_fdm
from compas.colors import Color

from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import cross_vectors
from compas.geometry import scale_vector
from compas_assembly.datastructures import Block
from compas_assembly.datastructures import Assembly


def ribs_from_pattern(form, width, depth):
    """Generate ribs with given width from a structural pattern and a floor depth.add()

    Parameters
    ----------
    form : FormDiagram
        The form diagram representing the structure
    width : float
        Width of the ribs in the pattern
    depth : float
        Depth of the floor (ribs will be extended until this height)

    Returns
    -------
    assembly: Assembly
        Asssembly with the ribs as blocks
    """

    assembly = Assembly()
    unitz = Vector(0, 0, 1)
    for edge in form.edges():
        u, v = edge
        length = form.edge_length(u, v)

        dir = form.edge_direction(u, v)
        vec = cross_vectors(dir, unitz)
        ptu = Point(*form.vertex_coordinates(u))
        ptv = Point(*form.vertex_coordinates(v))
        pt1, pt2 = ptu + scale_vector(vec, 1/2*width), ptu - scale_vector(vec, 1/2*width)
        pt3, pt4 = ptv + scale_vector(vec, 1/2*width), ptv - scale_vector(vec, 1/2*width)
        pt1_, pt2_ = pt1[:2] + [depth], pt2[:2] + [depth]
        pt3_, pt4_ = pt3[:2] + [depth], pt4[:2] + [depth]
        pts = [pt1, pt2, pt3, pt4, pt1_, pt2_, pt3_, pt4_]
        faces = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 2, 6, 4], [2, 3, 7, 6], [1, 3, 7, 5], [1, 0, 4, 5]]
        block = Block.from_vertices_and_faces(pts, faces)
        assembly.add_block(block)

    return assembly


xy_span = [[0.0, 8.1], [0.0, 8.1]]
xc, yc = sum(xy_span[0])/2, sum(xy_span[1])/2
floor_depth = 0.90
thickness = 0.05
density = 1.0
discretisation = 14
steps_loadpath = 4
width = 0.10

# form : FormDiagram = FormDiagram.create_cross_form(xy_span=xy_span, discretisation=discretisation)
form : FormDiagram = FormDiagram.create_fan_form(xy_span=xy_span, discretisation=discretisation)

# ------ COMPUTE LOADPATH OPTIMISATION IN A FEW LOOPS ------

for i in range(steps_loadpath):

    print('Compute loadpath iteration:', i)

    # Update tributary area of shell
    form.update_lumped_weights(thickness=thickness, density=density)

    # Add tributary area from ribs
    for edge in form.edges():
        u, v = edge
        length = form.edge_length(u, v)
        h1 = floor_depth - form.vertex_coordinates(u)[2]
        h2 = floor_depth - form.vertex_coordinates(v)[2]
        ph1 = -1 * h1 * width * 1/2 * length * density
        ph2 = -1 * h2 * width * 1/2 * length * density
        form.add_pz(u, ph1)
        form.add_pz(v, ph2)

    analysis = Analysis.create_lp_analysis(form, solver='CVXPY')
    analysis.set_up_optimiser()
    analysis.run()

    lp = form.loadpath()
    print('Loadpath of the solution before scallling is:', lp)

    vertical_equilibrium_fdm(form, zmax=floor_depth - thickness/2)
    lp = form.loadpath()
    print('Loadpath of the solution after scallling is:', lp)

form.overview_forces()

view = Viewer(form, show_grid=False)
view.draw_thrust()
view.show()

# Apply bounds
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    form.vertex_attribute(key, 'ub', floor_depth)
    form.vertex_attribute(key, 'lb', z - thickness/2)

assembly = ribs_from_pattern(form, width, floor_depth)

view = Viewer(form, show_grid=False)
view.settings['camera.target'] = [xc, yc, 0]
view.settings['camera.distance'] = 20
view.settings['opacity.shapes'] = 0.2
view.draw_thrust()
view.draw_shape()
view.draw_assembly(assembly, opacity=0.2)
# view.draw_mesh(mesh=upper_mesh, opacity=0.5, color=Color.from_rgb255(125, 125, 125))
view.show()

form_ad = '/Users/mricardo/compas_dev/me/floor/floor_form_' + form.parameters['type'] + '_depth_' + str(floor_depth) + '.json'
form.to_json(form_ad)

assembly_ad = '/Users/mricardo/compas_dev/me/floor/floor_asssembly_' + form.parameters['type'] + '_depth_' + str(floor_depth) + '.json'
assembly.to_json(assembly_ad)

print('Form Diagram create with the problem saved at:', form_ad)
print('Assembly create with the problem saved at:', assembly_ad)
