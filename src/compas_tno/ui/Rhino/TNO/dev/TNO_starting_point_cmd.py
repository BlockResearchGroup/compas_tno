from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino


__commandname__ = "TNO_starting_point"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']
    proxy = sc.sticky['TNO']['proxy']

    proxy.restart_server()

    objects = scene.find_by_name('Form')
    if not objects:
        compas_rhino.display_message("There is no FormDiagram in the scene.")
        return
    form = objects[0]

    objects = scene.find_by_name('Shape')
    if not objects:
        compas_rhino.display_message("There is no Shape in the scene.")
        return

    # add check to see if loads and

    proxy.package = 'compas_tno.problems'

    formdata = form.diagram.to_data()
    form.diagram.data, output = proxy.initialize_loadpath_proxy(formdata)

    # modify the basic visualisation after analysis
    form.settings['show.reactionvectors'] = True
    form.settings['show.cracks'] = True
    form.settings['show.forcepipes'] = True
    form.settings['show.reactionlabels'] = True

    message = output['status'] + ' - LP: ' + str(round(output['fopt'], 1)) + ' kN.m'

    compas_rhino.display_message(message)

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
