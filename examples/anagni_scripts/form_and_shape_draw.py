from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import Scene
import compas_rhino
import os

from compas_tno.rhino import SettingsForm

from compas_tno.rhino import FormObject
from compas_tno.rhino import ShapeObject
from compas_tno.rhino import OptimiserObject
from compas_tno.rhino import ForceObject


folder = '/Users/mricardo/compas_dev/me/anagni/'
objective = 'min'
title = 'top_vault_less_discr'
file_solution = os.path.join(folder, title +'_form_' + objective + '.json')
file_shape = os.path.join(folder, title + '_shape_' + objective + '.json')

form = FormDiagram.from_json(file_solution)
shape = Shape.from_json(file_shape)

scene = Scene()

guid_form = scene.add(form, name="Form", layer="TNO::FormDiagram")
guid_scene = scene.add(shape, name="Shape", layer="TNO::Shape")

object = scene.find(guid_form)

#object.settings['scale.forcepipes'] = 0.0003
object.settings['show.forcepipes'] = True
object.settings['show.reactionvectors'] = True
object.settings['show.reactionlabels'] = True
object.settings['show.cracks'] = True
#object.settings['scale.vectors'] = 0.001

scene.update()

#test = compas_rhino.rs.GetString("test", "test")

#SettingsForm.from_scene(scene, object_types=[FormObject, ShapeObject])

#scene.update()
