from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino

from compas_tno.rhino import SettingsForm
from compas_tno.rhino import FormObject
from compas_tno.rhino import ShapeObject
from compas_tno.rhino import OptimiserObject
from compas_tno.rhino import ForceObject


__commandname__ = "TNO_settings"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']
    if not scene:
        return

    SettingsForm.from_scene(scene, object_types=[FormObject, ForceObject, ShapeObject, OptimiserObject], global_settings=['TNO'])


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
