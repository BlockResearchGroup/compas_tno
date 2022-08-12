

from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
from compas_tno.shapes import Shape
from compas_tno.shapes import MeshDos
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer


objective = 'max'  # try 'max' 'Ecomp-linear'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']  # 't'
features = ['fixed']
starting_point = 'loadpath'

jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_rectangular/'
intra_file = jsonpath + '1-intrados-mesh.json'
extra_file = jsonpath + '1-extrados-mesh.json'

intrados = MeshDos.from_json(intra_file)
extrados = MeshDos.from_json(extra_file)

shape = Shape.from_meshes(intrados, extrados, data={'type': 'general', 'thk': 0.10, 't': 0.0})

# load meshes

for prob in ['A']:

    form_file = jsonpath + 'form-' + prob + '.json'

    form = FormDiagram.from_json(form_file)

    optimiser = Optimiser()
    optimiser.settings['library'] = solver
    optimiser.settings['solver'] = solver
    optimiser.settings['constraints'] = constraints
    optimiser.settings['variables'] = variables
    optimiser.settings['features'] = features
    optimiser.settings['objective'] = objective
    optimiser.settings['plot'] = False
    optimiser.settings['find_inds'] = True
    optimiser.settings['max_iter'] = 5000
    optimiser.settings['gradient'] = True
    optimiser.settings['jacobian'] = True
    optimiser.settings['printout'] = True
    optimiser.settings['starting_point'] = starting_point

    analysis = Analysis.from_elements(shape, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.set_up_optimiser()
    analysis.run()

    pzt = form.lumped_swt()
    thrust = form.thrust()
    print('Weight:', pzt)
    print('Thrust/Weight:', 100*round(thrust/pzt, 3))

    if optimiser.exitflag == 0:
        jsonsave = jsonpath + 'form-' + prob + '_' + objective +'.json'
        form.to_json(jsonsave)

    viewer = Viewer(form, shape)
    viewer.settings['camera.show.grid'] = False
    viewer.settings['camera.distance'] = 10
    viewer.settings['camera.target'] = [1.5, 1, 0]
    viewer.settings['scale.reactions'] = 0.1
    viewer.settings['force.anchor'] = [5, 5, 0]
    viewer.settings['force.scale'] = 0.1
    viewer.draw_shape()
    viewer.show()
    # viewer.show_solution()

    plotter = TNOPlotter(form)
    plotter.show_solution()
