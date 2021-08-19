from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import Scene

data = {'type': 'cross_fd', 'xy_span': [[0.0, 10.0], [0.0, 10.0]], 'discretisation': 10, 'fix': 'corners'}

form = FormDiagram.from_library(data)

scene = Scene()

scene.add(form, name='Form', layer='TNO::FormDiagram')
scene.update()

print('ok')
