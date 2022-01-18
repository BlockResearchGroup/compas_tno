from compas_tno.diagrams import FormDiagram
import compas_tno
# from compas_tno.viewers import animation_from_optimisation
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.viewers import Viewer

# ---- MAKE THE VIDEO ----

DATA_FORM = compas_tno.get('form.json')
# DATA_XFORM = compas_tno.get('Xform.json')

form = FormDiagram.from_json(DATA_FORM)
# animation_from_optimisation(form, DATA_XFORM, interval=150)

force = reciprocal_from_form(form)

from compas.geometry import Translation
from compas.geometry import Scale

S = Scale.from_factors([1/20, 1/20, 1/20])
T = Translation.from_vector([0.0, -10, 0])

force = force.transformed(S)
force = force.transformed(T)

print(force)
print(type(force))

viewer = Viewer(form)
viewer.view_force(force)
viewer.show_solution()
