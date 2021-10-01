from compas_tno.shapes import Shape
from compas_tno.rhino import Scene
import compas_rhino

type_structure = 'dome'
discretisation = [20, 16]

# type_structure = 'crossvault'
# discretisation = 20

# type_structure = 'pavillionvault'
# discretisation = 20

#type_structure = 'pointed_crossvault'
#discretisation = 20
#hc = 6.5

data_shape = {'xy_span': [[0.0, 6.0], [0.0, 10.0]], 'center': [5.0, 5.0], 'radius': 5.0, 'type': type_structure, 'thk': 0.5, 'discretisation': discretisation, 'hc': hc, 't': 0.0}
shape = Shape.from_library_proxy(data_shape)

scene = Scene()

scene.add(shape, name="Shape", layer="TNO::Shape")
scene.update()