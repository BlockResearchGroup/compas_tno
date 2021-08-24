from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino
from compas_rhino.geometry import RhinoMesh
from compas_tno.shapes import Shape


__commandname__ = "TNO_shape_from_meshes"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']

    answer = compas_rhino.rs.GetString("Select a method to input the shape. From", "Cancel", ["IntraExtrados", "IntraExtradosMiddle", "Middle", "Cancel"])
    if not answer:
        return
    if answer == "Cancel":
        return

    middle = None

    warning = 'Warning: Select only one mesh.'

    if answer == "IntraExtrados" or answer == "IntraExtradosMiddle":
        selection = compas_rhino.select_meshes(message='Select Intrados Mesh')
        if not selection:
            return
        if len(selection) > 1:
            compas_rhino.display_message(warning)
        guid_intrados = selection[0]
        intrados_RM = RhinoMesh.from_guid(guid_intrados)
        intrados = intrados_RM.to_compas()

        compas_rhino.rs.UnselectObjects(selection)

        selection = compas_rhino.select_mesh(message='Select Extrados Mesh')
        if not selection:
            return
        if len(selection) > 1:
            compas_rhino.display_message(warning)
        guid_extrados = selection[0]
        extrados_RM = RhinoMesh.from_guid(guid_extrados)
        extrados = extrados_RM.to_compas()

        compas_rhino.rs.UnselectObjects(selection)

    if answer == "Middle" or answer == "IntraExtradosMiddle":
        selection = compas_rhino.select_mesh(message='Select Middle Mesh')
        if not selection:
            return
        if len(selection) > 1:
            compas_rhino.display_message(warning)
        guid_middle = selection[0]
        middle_RM = RhinoMesh.from_guid(guid_middle)
        middle = middle_RM.to_compas()

    thk = compas_rhino.rs.GetReal("Input thickness value (estimate)", 0.50)
    if not thk:
        return

    shape = Shape.from_meshes(intrados, extrados, middle, data={'type': 'general', 't': 0.0, 'thk': thk})

    scene.add(shape, name='Shape', layer='TNO::Shape')
    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
