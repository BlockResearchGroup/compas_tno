from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino


__commandname__ = "TNO_run_optimisation"


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

    objects = scene.find_by_name('Optimiser')
    if not objects:
        compas_rhino.display_message("Select First the optimisation settings.")
        return
    optimiser = objects[0]

    # add check to see if loads are applied and if bounds on q are set

    proxy.package = 'compas_tno.problems'

    formdata = form.diagram.to_data()
    shapedata = shape.shape.to_data()
    optimiser.optimiser.message = None
    optimiserdata = optimiser.optimiser.to_data()

    shapedata, formdata, optimiserdata = proxy.run_NLP_proxy(shapedata, formdata, optimiserdata)

    form.diagram.data = formdata
    shape.shape.data = shapedata
    optimiser.optimiser.data = optimiserdata

    message = optimiser.optimiser.message + ' fopt: ' + str(round(optimiser.optimiser.fopt, 2))
    compas_rhino.display_message(message)

    # modify the basic visualisation after analysis
    form.settings['show.reactionvectors'] = True
    form.settings['show.cracks'] = True
    form.settings['show.forcepipes'] = True
    form.settings['show.reactionlabels'] = True

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
