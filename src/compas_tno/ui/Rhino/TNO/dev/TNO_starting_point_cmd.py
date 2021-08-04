from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.diagrams import FormGraph
from compas_tno.diagrams import FormDiagram


__commandname__ = "TNO_starting_point"


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


    # scene = sc.sticky['TNO']['scene']

    # guids = compas_rhino.select_lines(message='Select Form Diagram Lines')
    # if not guids:
    #     return

    # lines = compas_rhino.get_line_coordinates(guids)
    # graph = FormGraph.from_lines(lines)

    # if not graph.is_planar_embedding():
    #     compas_rhino.display_message('The graph is not planar. Therefore, a form diagram cannot be created.')
    #     return

    # form = FormDiagram.from_graph(graph)

    # scene.purge()
    # scene.add(form, name='Form', layer='TNO::FormDiagram')
    # scene.update()
    # scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
