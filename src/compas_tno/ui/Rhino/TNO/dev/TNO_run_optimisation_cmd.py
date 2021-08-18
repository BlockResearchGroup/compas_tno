from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.diagrams import FormGraph
from compas_tno.diagrams import FormDiagram


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
    shapedata = shape.shape.datashape  # WIP
    optimiserdata = optimiser.optimiser.settings

    optimiserdata['solver'] = 'SLSQP'

    shapedata, formdata, optimiserdata = proxy.run_NLP_proxy(shapedata, formdata, optimiserdata)

    form.diagram.data = formdata

    print(optimiserdata)

    message = optimiserdata['status'] + ' fopt: ' + str(round(optimiserdata['fopt'], 2))
    compas_rhino.display_message(message)

    form.settings['show.cracks'] = True

    scene.update()
    scene.save()



# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
