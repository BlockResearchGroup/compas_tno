from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import scriptcontext as sc

import compas_rhino
from compas_tno.optimisers import Optimiser


__commandname__ = "TNO_optimization_form"


def RunCommand(is_interactive):

    if 'TNO' not in sc.sticky:
        compas_rhino.display_message('TNO has not been initialised yet.')
        return

    scene = sc.sticky['TNO']['scene']

    data = {}

    data['printout'] = True

    # objective
    obj = compas_rhino.rs.GetString("Optimisation Objective", "min", ["min", "max", "t", "Cancel"])
    if not obj or obj == "Cancel":
        return
    data['objective'] = obj

    # variables
    options = ['q', 'zb', 'xyb', 't']
    items = ("ForceDensitiesQ", "False", "True"), ("HeightSupportsZb", "False", "True"), ("PlanarSupportsXyb", "False", "True"), ("Thickness", "False", "True")
    results = compas_rhino.rs.GetBoolean("Select Variables", items, (True, True, False, False))
    if not results:
        return
    data['variables'] = [i for indx, i in enumerate(options) if results[indx]]

    # constraints
    options = ['funicular', 'envelope', 'reac_bounds', 'envelopexy']
    items = ("Funicular", "False", "True"), ("Envelope", "False", "True"), ("ReactionBounds", "False", "True"), ("EnvelopeXy", "False", "True")
    results = compas_rhino.rs.GetBoolean("Select Constraints", items, (True, True, False, False))
    if not results:
        return
    data['constraints'] = [i for indx, i in enumerate(options) if results[indx]]

    # features
    options = ['fixed', 'sym']
    items = ("FixedProjection", "False", "True"), ("Symmetry", "False", "True")
    results = compas_rhino.rs.GetBoolean("Select Features", items, (True, False))
    if not results:
        return
    data['features'] = [i for indx, i in enumerate(options) if results[indx]]

    # solver
    lib = compas_rhino.rs.GetString("Select Library", "SLSQP", ["SLSQP", "IPOPT", "Cancel"])
    if not lib or lib == "Cancel":
        return
    data['library'] = lib
    data['solver'] = lib

    # derivatives
    deriv = compas_rhino.rs.GetString("Use analytical derivatives", "True", ["True", "False", "Cancel"])
    if deriv == "True":
        data['gradient'] = True
        data['jacobian'] = True

    # max iter
    max_iter = compas_rhino.rs.GetInteger("Maximum of Iterations", 500)
    if not max_iter:
        return
    data['max_iter'] = max_iter

    # starting point
    sp = compas_rhino.rs.GetString("Select Starting point", "current", ["current", "loadpath", "ParalelliseTna", "Cancel"])
    if not sp or sp == "Cancel":
        return
    data['starting_point'] = sp

    # qmax qmin
    qmin = compas_rhino.rs.GetReal("Assign qmin", -1e+4)
    if not qmin:
        pass
    else:
        data['qmin'] = float(qmin)
    qmax = compas_rhino.rs.GetReal("Assign qmax", +1e-8)
    if not qmax:
        pass
    else:
        data['qmax'] = float(qmax)

    objects = scene.find_by_name('Optimiser')
    if not objects:
        optimiser = Optimiser()
        optimiser.settings = data
        scene.add(optimiser, name='Optimiser', layer=None)
        objects = scene.find_by_name('Optimiser')
        optimiserobject = objects[0]
    else:
        optimiserobject = objects[0]
        optimiser = optimiserobject.optimiser
        optimiser.settings = data

    optimiserobject.update_object_from_optimiser()

    scene.update()
    scene.save()


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':

    RunCommand(True)
