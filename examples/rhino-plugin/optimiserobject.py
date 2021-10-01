import compas_rhino
from compas_tno.optimisers import Optimiser
from compas_tno.rhino import OptimiserObject
from compas_tno.rhino import SettingsForm
from compas_tno.rhino import Scene

scene = Scene()

datadefault = {'printout': True, 'features': ['blabla'], 'constraints': ['funicular', 'envelope'], 'gradient': True, 'objective': 'max', 'qmax': 1e-08, 'variables': ['q', 'zb'], 'max_iter': 500, 'solver': 'SLSQP', 'jacobian': True, 'library': 'SLSQP', 'starting_point': 'current', 'qmin': -10000.0}

optimiser = Optimiser()
optimiser.settings = datadefault

scene.add(optimiser, name='Optimiser', layer=None)
objects = scene.find_by_name('Optimiser')
optimiserobject = objects[0]
optimiserobject.update_object_from_optimiser()

SettingsForm.from_scene(scene, object_types=[OptimiserObject])

print(optimiserobject.optimiser.settings)