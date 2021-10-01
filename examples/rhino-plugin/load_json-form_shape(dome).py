import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.rhino import Scene
from compas_tno.rhino import SettingsForm

from compas_tno.rhino import FormObject
from compas_tno.rhino import ShapeObject
from compas_tno.rhino import OptimiserObject
from compas_tno.rhino import ForceObject


json_shape = '/Users/mricardo/compas_dev/compas_tno/data/shape-dome.json'
json_form = '/Users/mricardo/compas_dev/compas_tno/data/form-dome.json'

form = FormDiagram.from_json(json_form)
shape = Shape.from_json(json_shape)

SETTINGS = {
    'TNO': {
        'autoupdate': True,
    }
}

scene = Scene(settings=SETTINGS)

scene.add(form, name="Form", layer="TNO::Form")
scene.add(shape, name="Shape", layer="TNO::Shape")

objects = scene.find_by_name("Form")

formobject = objects[0]
formobject.settings['scale.forcepipes'] = 0.0005
formobject.settings['show.edges'] = False
formobject.settings['show.reactionvectors'] = True
formobject.settings['show.cracks'] = True
formobject.settings['show.forcepipes'] = True
formobject.settings['show.reactionlabels'] = True

scene.update()

#SettingsForm.from_scene(scene, object_types=[FormObject, ShapeObject, OptimiserObject, ForceObject], global_settings=['TNO'])

#scene.update()

