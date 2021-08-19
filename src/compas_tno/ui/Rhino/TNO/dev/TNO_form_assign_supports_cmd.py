from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

__commandname__ = "TNO_form_assign_supports"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']

    objects = scene.find_by_name('Form')
    if not objects:
        compas_rhino.display_message("There is no FormDiagram in the scene.")
        return
    form = objects[0]

    form.diagram.vertices_attribute('is_fixed', False)

    show_vertices = form.settings['show.vertices']
    form.settings['show.vertices'] = True if not form.settings['show.vertices'] else show_vertices

    scene.update()

    vertices = form.select_vertices("Fix selected vertices (unfix all others)")
    if not vertices:
        return

    form.diagram.vertices_attribute('is_fixed', True, keys=vertices)
    form.settings['show.vertices'] = show_vertices

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
