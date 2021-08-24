import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.rhino import Scene

json_shape = '/Users/mricardo/compas_dev/compas_tno/data/shape.json'
json_form = '/Users/mricardo/compas_dev/compas_tno/data/form.json'

form = FormDiagram.from_json(json_form)
shape = Shape.from_json(json_shape)

scene = Scene()

scene.add(form, name="Form", layer="TNO::Form")
scene.add(shape, name="Shape", layer="TNO::Shape")

objects = scene.find_by_name("Form")

formobject = objects[0]
formobject.settings['show.reactionvectors'] = True
formobject.settings['show.loadvectors'] = True
formobject.settings['show.cracks'] = True
formobject.settings['show.forcepipes'] = True
formobject.settings['show.reactionlabels'] = True


scene.update()