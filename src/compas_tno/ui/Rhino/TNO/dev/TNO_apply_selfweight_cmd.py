from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape


__commandname__ = "TNO_apply_selfweight"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']
    proxy = sc.sticky['TNO']['proxy']

    objects = scene.find_by_name('Form')
    if not objects:
        compas_rhino.display_message("There is no FormDiagram in the scene.")
        return
    form = objects[0]

    objects = scene.find_by_name('Shape')
    if not objects:
        compas_rhino.display_message("There is no Shape in the scene.")
        return
    shape = objects[0]

    proxy.package = 'compas_tno.utilities'

    formdata = form.diagram.to_data()
    shapedata = shape.shape.datashape  # WIP
    form.diagram.data = proxy.apply_selfweight_from_shape_proxy(formdata, shapedata)

    # form = FormDiagram.from_data(formdata)

    # scene.purge()
    # scene.add(form, name='Form', layer='TNO::FormDiagram')
    # scene.add(shape, name='Shape', layer='TNO::Shape')

    form.settings['show.vertexloads'] = True

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)