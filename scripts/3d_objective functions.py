from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
import compas_tno
from compas.geometry import Line
from compas_view2.shapes import Arrow
from compas.geometry import normalize_vector
from compas.geometry import norm_vector
from numpy import array
from numpy import zeros
# from compas_tno.utilities import apply_envelope_from_shape

delta = 1.0
span = 10.0
xspan = yspan = [0.0, span]
xspan_vault = yspan_vault = [- delta, span + delta]
thk = 0.50

data_vault = {
    'type': 'crossvault',
    'xy_span': [xspan_vault, yspan_vault],
    'thk': thk,
    'discretisation': 10,
    't': 0.0
}

data_diagram = {
    'type': 'fan_fd',
    'xy_span': [xspan, yspan],
    'thk': thk,
    'discretisation': 10,
    'fix': 'corners'
}

vault = Shape.from_library(data_vault)
form = FormDiagram.from_library(data_diagram)

from compas_plotters import Plotter
plotter = Plotter()
plotter.fontsize = 6
artist = plotter.add(form)
artist.draw_vertexlabels()
plotter.zoom_extents()
plotter.show()

objective = 'bestfit'  # try 'max' 'Ecomp-linear'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']  # 'lambdv', 't'
features = ['fixed', 'sym']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['objective'] = objective
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = False
optimiser.settings['max_iter'] = 500
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['printout'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['sym_loads'] = False

if objective == 'max_load':
    # loaded_node = 58
    loaded_node = 92
    max_load_mult = 2000.0
    n = form.number_of_vertices()
    pzv = zeros((n, 1))
    pzv[loaded_node] = -1.0

    optimiser.settings['max_lambd'] = max_load_mult
    optimiser.settings['load_direction'] = pzv

if objective == 'Ecomp-linear':
    Xc = [5.0, 5.0, 0.0]
    vector_supports = []
    sign = +1  # +1 for outwards / -1 for inwards

    for key in form.vertices_where({'is_fixed': True}):
        x, y, z = form.vertex_coordinates(key)
        xi, yi, zi = vault.intrados.vertex_coordinates(key)

        # dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), sign*(z - Xc[2])])  # 4 corners

        if x > Xc[0] and y < Xc[1]:             # vertical settlement
            dXbi = normalize_vector([0, 0, -1])
        else:
            dXbi = [0, 0, 0]

        vector_supports.append(dXbi)
        form.vertex_attribute(key, 'dXb', dXbi)

    dXb = array(vector_supports)
    print(dXb)

    optimiser.settings['support_displacement'] = dXb

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()

vault0 = Shape.from_formdiagram_and_attributes(form)

analysis.run()

path = compas_tno.get('')
form_path = path + '/form-' + objective + '.json'
analysis_path = path + '/analysis-' + objective + '.json'
form.to_json(form_path)
print('Form saved to:', form_path)
# analysis.to_json(analysis_path)

vault2 = Shape.from_formdiagram_and_attributes(form)
viewer = Viewer(form, vault2)
viewer.settings['camera.show.grid'] = False
viewer.settings['camera.distance'] = 35

viewer.draw_thrust()
# viewer.draw_middle_shape()

if objective == 't':
    viewer.settings['color.mesh.intrados'] = (255, 150, 150)
    viewer.settings['color.mesh.extrados'] = (255, 150, 150)
    viewer.draw_mesh(vault0.intrados, show_edges=False, color=(125, 125, 125))
    viewer.draw_mesh(vault0.extrados, show_edges=False, color=(125, 125, 125))
else:
    viewer.draw_cracks()

viewer.draw_shape()

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
