from compas_tno.shapes import Shape
from compas_tno.rhino import Scene
from compas_tno.diagrams import FormDiagram
from compas.rpc import Proxy

data_shape = {'xy_span': [[0.0, 10.0], [0.0, 10.0]], 'type': 'crossvault', 'thk': 0.5, 'discretisation': 10, 't': 0.0}
shape = Shape.from_library_proxy(data_shape)

data_form = {'type': 'cross_fd', 'xy_span': [[0.0, 10.0], [0.0, 10.0]], 'discretisation': 10, 'fix': 'corners'}
form = FormDiagram.from_library(data_form)

scene = Scene()
proxy = Proxy()

proxy.restart_server()

proxy.package = 'compas_tno.utilities'

formdata = form.to_data()
shapedata = data_shape
formdata = proxy.apply_selfweight_from_shape_proxy(formdata, shapedata)
formdata = proxy.apply_envelope_from_shape_proxy(formdata, shapedata)
formdata = proxy.apply_bounds_on_q_proxy(formdata)

form = FormDiagram.from_data(formdata)

scene.add(shape, name="Shape", layer="TNO::Shape")
scene.add(form, name='Form', layer='TNO::FormDiagram')
scene.update()

print('Part I - ok')

# proxy = Proxy()

objects = scene.find_by_name('Form')
formobject = objects[0]

objects = scene.find_by_name('Shape')
shapeobject = objects[0]

proxy.package = 'compas_tno.problems'

formdata = form.to_data()
formobject.diagram.data = proxy.initialize_loadpath_proxy(formdata)

scene.update()

print('Part II - ok')