from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.shapes import Shape


__commandname__ = "TNO_shape_from_library"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']

    answer = compas_rhino.rs.GetString("Base geometry of the Shape", "Cancel", ["Circular", "Rectangular", "Cancel"])
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

        layout = compas_rhino.rs.GetString("Set the type of shape (See docs)", "Cancel", ["dome", "dome_polar", "dome_spr", "Cancel"])
        if not layout:
            return
        if layout == "Cancel":
            return
        data['type'] = layout

        thk = compas_rhino.rs.GetReal("Select the thickness", 0.5)
        if not thk:
            return
        data['thk'] = thk

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

        t = compas_rhino.rs.GetReal("Define t (See docs.)", 0.0)
        if t:
            data['t'] = t

    if answer == "Rectangular":
        rect = compas_rhino.rs.GetRectangle(mode=1, prompt1='Select the first corner of the diagram', prompt2='Select the second corner of the diagram')
        if not rect:
            return
        xy_coords = [[], []]
        for corner in rect:
            xy_coords[0].append(corner[0])
            xy_coords[1].append(corner[1])

        data['xy_span'] = [[min(xy_coords[0]), max(xy_coords[0])], [min(xy_coords[1]), max(xy_coords[1])]]

        layout = compas_rhino.rs.GetString("Set the type of shape (See docs)", "Cancel", [
                                           "crossvault", "pavillionvault", "parabolic_shell", "pointed_crossvault", "domicalvault", "Cancel"])
        if not layout:
            return
        if layout == "Cancel":
            return
        data['type'] = layout

        thk = compas_rhino.rs.GetReal("Select the thickness", 0.5)
        if not thk:
            return
        data['thk'] = thk

        discr = compas_rhino.rs.GetInteger(message="Select the discretisation", number=10)
        if not discr:
            return
        data['discretisation'] = discr

        t = compas_rhino.rs.GetReal("Select parameter t", 0.0)
        data['t'] = t

        if data['type'] == 'pointed_crossvault':
            min_hc = max((data['xy_span'][0][1] - data['xy_span'][0][0])/2, (data['xy_span'][1][1] - data['xy_span'][1][0])/2)
            hc = compas_rhino.rs.GetReal("Height at midspan (>{0})".format(min_hc), min_hc)
            if not hc:
                return
            data['hc'] = float(max(hc, min_hc))

            add_h = compas_rhino.rs.GetString("Set heights on edge and rubble?", "No", ["Yes", "No", "Cancel"])
            if add_h == "Yes":
                he = compas_rhino.rs.GetReal("Height at edge (<{0})".format(min_hc), min_hc)
                if he:
                    data['he'] = he
                hm = compas_rhino.rs.GetReal("Height at rubble (<{0})".format(min_hc), min_hc)
                if hm:
                    data['hm'] = hm
            else:
                data['hm'] = None
                data['he'] = None

    scene.clear()

    shape = Shape.from_library_proxy(data)

    objects = scene.find_by_name('Shape')
    if not objects:
        scene.add(shape, name='Shape', layer='TNO::Shape')
    else:
        shapeobject = objects[0]
        shapeobject.shape = shape

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
