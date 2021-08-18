
from compas_tno.rhino import Scene
from compas_tno.optimisers import Optimiser


scene = Scene()
data = {'constraints': ['funicular', 'envelope'], 'objective': 'min', 'qmax': 1e-08, 'variables': ['q', 'zb'], 'max_iter': 500, 'fetures': ['fixed'], 'library': 'SLSQP', 'starting_point': 'loadpath', 'qmin': -10000.0}

optimiser = Optimiser()
optimiser.settings = data

scene.add(optimiser, name='Optimiser', layer=None)
scene.update()
