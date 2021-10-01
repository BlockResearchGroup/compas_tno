from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import Scene
import compas_rhino

from compas_tno.rhino import SettingsForm

from compas_tno.rhino import FormObject
from compas_tno.rhino import ShapeObject
from compas_tno.rhino import OptimiserObject
from compas_tno.rhino import ForceObject


type_structure = 'dome'
discretisation = [20, 16]

#data_shape = {'center': [5.0, 5.0], 'radius': 5.0, 'type': type_structure, 'thk': 0.5, 'discretisation': discretisation, 't': 0.0}
#shape = Shape.from_library_proxy(data_shape)

jsonp = '/Users/mricardo/compas_dev/compas_tno/data/form-dome.json'
form = FormDiagram.from_json(jsonp)

scene = Scene()

guid = scene.add(form, name="Form", layer="TNO::FormDiagram")
# scene.add(shape, name="Shape", layer="TNO::Shape")

scene.update()

object = scene.find(guid)

object.settings['scale.forcepipes'] = 0.0003
object.settings['show.forcepipes'] = True
object.settings['show.reactionvectors'] = True
object.settings['scale.vectors'] = 0.001

test = compas_rhino.rs.GetString("test", "test")

SettingsForm.from_scene(scene, object_types=[FormObject])

scene.update()
