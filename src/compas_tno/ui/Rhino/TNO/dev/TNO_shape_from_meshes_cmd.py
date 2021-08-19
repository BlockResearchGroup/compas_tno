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

    if answer == "IntraExtrados" or answer == "IntraExtradosMiddle":
        guid_intrados = compas_rhino.select_mesh(message='Select Intrados Mesh')
        if not guid_intrados:
            return
        intrados_RM = RhinoMesh.from_guid(guid_intrados)
        intrados = intrados_RM.to_compas()

        guid_extrados = compas_rhino.select_mesh(message='Select Intrados Mesh')
        if not guid_extrados:
            return
        extrados_RM = RhinoMesh.from_guid(guid_extrados)
        extrados = extrados_RM.to_compas()

    if answer == "Middle" or answer == "IntraExtradosMiddle":
        guid_middle = compas_rhino.select_mesh(message='Select Intrados Mesh')
        if not guid_middle:
            return
        middle_RM = RhinoMesh.from_guid(guid_middle)
        middle = middle_RM.to_compas()

    shape = Shape.from_meshes(intrados, extrados, middle, data=None)

    scene.add(shape, name='Shape', layer='TNO::Shape')
    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
