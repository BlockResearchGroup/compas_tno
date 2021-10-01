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

min = 'radial_spaced_discr_8_20_min_t=50.json'
max = 'radial_spaced_discr_8_20_max_t=50.json'
min_thk = 'radial_spaced_discr_8_20_min_t=207.json'

max_flower = 'flower_discr_8_20_max_t=50.json'
min_flower = 'flower_discr_8_20_min_t=50.json'
min_thk_flower = 'flower_discr_8_20_min_t=365.json'

#end = min_thk
#file_solution = '/Users/mricardo/compas_dev/me/minmax/dome/radial_spaced/'+ end

end = min_thk_flower
file_solution = '/Users/mricardo/compas_dev/me/minmax/dome/flower/'+ end

# vaults

end = 'cross_fd_discr_20_min_t=50.json'
#end = 'cross_fd_discr_20_max_t=50.json'
#end = 'cross_fd_discr_20_min_t=304.json'
file_solution = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/cross_fd/' + end

end = form.json
file_solution = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/cross_fd/' + end


#end = 'fan_fd_discr_16_min_t=50.json'
#end = 'fan_fd_discr_16_max_t=50.json'
#end = 'fan_fd_discr_16_min_t=42.json'
#file_solution = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/fan_fd/' + end


form = FormDiagram.from_json(file_solution)

shape = Shape.from_formdiagram_and_attributes(form)

for edge in form.edges():
    qi = form.edge_attribute(edge, 'q')
    form.edge_attribute(edge, 'q', -1 * qi)

#for edge in form.edges_on_boundary():
#    form.edge_attribute(edge, '_is_edge', False)

scene = Scene()

guid_form = scene.add(form, name="Form", layer=end)

guid_shape = scene.add(shape, name="Shape", layer="TNO::Shape")

object = scene.find(guid_form)

#object.settings['scale.forcepipes'] = 0.0015
object.settings['scale.forcepipes'] = 0.01
object.settings['show.supports'] = False
object.settings['show.forcepipes'] = True
object.settings['show.reactionvectors'] = False
object.settings['show.reactionlabels'] = False
object.settings['show.cracks'] = True
#object.settings['scale.vectors'] = 0.001

scene.update()

#test = compas_rhino.rs.GetString("test", "test")

#SettingsForm.from_scene(scene, object_types=[FormObject, ShapeObject])

#scene.update()
