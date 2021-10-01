from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.diagrams import FormDiagram


__commandname__ = "TNO_form_from_library"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']

    answer = compas_rhino.rs.GetString("Geometry of the Form Diagram", "Cancel", ["Circular", "Rectangular", "Cancel"])
    if not answer:
        return
    if answer == "Cancel":
        return

    data = {}
    if answer == "Circular":
        center = compas_rhino.rs.GetPoint(message='Select the center of the diagram')
        if not center:
            return
        data['center'] = [center[0], center[1]]

        radius = compas_rhino.rs.GetReal("Set Radius", 5.0)
        if not radius:
            radius = compas_rhino.rs.GetDistance(first_pt=center, second_pt_msg='Select Second Point to find radius')
        if not radius:
            return
        data['radius'] = radius

        layout = compas_rhino.rs.GetString("Set the layout of the diagram (See docs)", "Cancel", ["radial_fd", "radial_spaced_fd", "spiral_fd", "Cancel"])
        if not layout:
            return
        if layout == "Cancel":
            return
        data['type'] = layout

        discr_1 = compas_rhino.rs.GetInteger(message="Radial discretisaton (#hoops)", number=8)
        if not discr_1:
            return
        discr_2 = compas_rhino.rs.GetInteger(message="Meridian discretisaton", number=12)
        if not discr_2:
            return
        data['discretisation'] = [discr_1, discr_2]

        oculus = compas_rhino.rs.GetString("Add Oculus?", "False", ["True", "False", "Cancel"])
        if oculus == "True":
            radius_oc = compas_rhino.rs.GetReal("Radius of oculus (<{0})".format(radius), radius/4)
            if radius_oc:
                data['r_oculus'] = radius_oc

        if layout in ["radial_fd", "radial_spaced_fd"]:
            diagonals = compas_rhino.rs.GetString("Add Diagonals?", "False", ["True", "False", "Cancel"])
            if diagonals == "True":
                data['diagonal'] = True
                diagonals_direction = compas_rhino.rs.GetString("Choose direction", "right", ["right", "left"])
                data['partial_diagonal'] = diagonals_direction

    if answer == "Rectangular":
        rect = compas_rhino.rs.GetRectangle(mode=1, prompt1='Select the first corner of the diagram', prompt2='Select the second corner of the diagram')
        if not rect:
            return
        xy_coords = [[], []]
        for corner in rect:
            xy_coords[0].append(corner[0])
            xy_coords[1].append(corner[1])

        data['xy_span'] = [[min(xy_coords[0]), max(xy_coords[0])], [min(xy_coords[1]), max(xy_coords[1])]]

        layout = compas_rhino.rs.GetString("Set the layout of the diagram (See docs)", "Cancel", ["cross_fd", "fan_fd", "cross_with_diagonal", "ortho", "Cancel"])
        if not layout:
            return
        if layout == "Cancel":
            return
        data['type'] = layout

        discr = compas_rhino.rs.GetInteger(message="Select the discretisation", number=10)
        if not discr:
            return
        data['discretisation'] = discr

        fix = compas_rhino.rs.GetString("Select fixity", "Cancel", ["all", "corners", "Cancel"])
        if not fix:
            return
        if fix == "Cancel":
            return
        data['fix'] = fix

    scene.clear()

    form = FormDiagram.from_library(data)

    objects = scene.find_by_name('Form')
    if not objects:
        scene.add(form, name='Form', layer='TNO::FormDiagram')
    else:
        formobject = objects[0]
        formobject.diagram = form

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
