from cvxpy import vec
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
from compas_assembly.datastructures import Assembly

from compas_view2.objects import Arrow

from numpy import array
from compas.geometry import normalize_vector
from compas.geometry import norm_vector


xy_span = [[0.0, 8.1], [0.0, 8.1]]
xc, yc = sum(xy_span[0])/2, sum(xy_span[1])/2
floor_depth = 0.90
thickness = 0.05
density = 1.0
discretisation = 14
steps_loadpath = 4
width = 0.10
solver = 'IPOPT'

type_fd = 'cross_fd'

form_ad = '/Users/mricardo/compas_dev/me/floor/floor_form_' + type_fd + '_depth_' + str(floor_depth) + '_optimisation_min.json'
assembly_ad = '/Users/mricardo/compas_dev/me/floor/floor_assembly_' + type_fd + '_depth_' + str(floor_depth) + '.json'

print('Form Diagram loaded from the problem saved at:', form_ad)
print('Assembly loaded from the problem saved at:', assembly_ad)

form = FormDiagram.from_json(form_ad)
assembly = Assembly.from_json(assembly_ad)

vector_supports = []
bases = []
for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    bases.append([x, y, z])
    dXbi = [0, 0, 0]
    if x > xc and y > yc:
        dXbi = normalize_vector([0, 0, -1])  # 4 corners
    vector_supports.append(dXbi)
dXb = array(vector_supports)

print('Initial LP:', form.loadpath())

shape = Shape.from_formdiagram_and_attributes(form)
analysis = Analysis.create_compl_energy_analysis(form, shape, solver=solver, starting_point='current', printout=True, support_displacement=dXb)
analysis.set_up_optimiser()
analysis.run()

thrust, swt = form.thrust(), form.lumped_swt()
print('Thrust/Weight:', round(thrust/swt, 3))

view = Viewer(form, show_grid=False)
view.settings['camera.target'] = [xc, yc, 0]
view.settings['camera.distance'] = 20
view.settings['opacity.shapes'] = 0.2

for i in range(len(vector_supports)):
    if norm_vector(vector_supports[i]) > 0.001:
        arrow = Arrow(bases[i], vector_supports[i])
        view.add(arrow, color=(0, 0, 0))

view.draw_thrust()
view.draw_cracks()
view.draw_shape()
view.draw_assembly(assembly, opacity=0.2)
# view.draw_mesh(mesh=upper_mesh, opacity=0.5, color=Color.from_rgb255(125, 125, 125))
view.show()
