from compas_tno.shapes import Shape
from compas_tno.rhino import Scene
from compas.rpc import Proxy
import compas_rhino

data_shape = {'xy_span': [[0.0, 10.0], [0.0, 10.0]], 'type': 'crossvault', 'thk': 0.5, 'discretisation': 10, 't': 0.0}
shape = Shape.from_library_proxy(data_shape)

scene = Scene()
proxy = Proxy()

scene.add(shape, name="Shape", layer="TNO::Shape")
scene.update()

objects = scene.find_by_name('Shape')
shapeobject = objects[0]

#print(shapeobject.intrados)
#print(shapeobject.shape)
#print(shapeobject.guids)

resp = compas_rhino.rs.GetString("Continue", "True", ["True", "Cancel"])

data_shape['thk'] = 3.00

shape2 = Shape.from_library_proxy(data_shape)
shapeobject.shape = shape2

scene.update()

