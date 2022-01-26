from compas_tno.diagrams import FormDiagram
import compas_tno
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.shapes.shape import Shape
from compas_tno.viewers import animation_from_optimisation
from compas_tno.viewers import animation_from_section
from compas_tno.viewers import Viewer

# ---- MAKE THE VIDEO ----

DATA_FORM = compas_tno.get('form.json')
DATA_SHAPE = compas_tno.get('shape.json')
DATA_XFORM = compas_tno.get('Xform.json')
DATA_XFORCE = compas_tno.get('Xforce.json')

form = FormDiagram.from_json(DATA_FORM)
shape = Shape.from_library(
    {
        'type': 'dome_polar',
        'thk': 0.50,
        # 'thk': 0.205,
        'discretisation': [20, 40],
        't': 0.0,
        'center': [5.0, 5.0],
        'radius': 5.0,
    })

force = reciprocal_from_form(form)

viewer = Viewer(form, shape=shape)
viewer.view_thrust()
viewer.view_force()
viewer.view_shape()
viewer.show()

SETTINGS = viewer.settings

# settings dome
SETTINGS['camera.show.grid'] = False
SETTINGS['camera.show.axis'] = False
SETTINGS['camera.distance'] = 20
SETTINGS['camera.target'] = [0, 0, 5]
SETTINGS['force.anchor'] = [10, -6, 0]
SETTINGS['force.scale'] = 0.04

# animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, settings=SETTINGS, record=True, interval=150)

animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, shape=shape, settings=SETTINGS, record=True, interval=20)

# animation_from_section(form, DATA_XFORM)
