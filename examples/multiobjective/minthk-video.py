from compas_tno.diagrams import FormDiagram
import compas_tno
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.viewers import animation_from_optimisation

# ---- MAKE THE VIDEO ----

DATA_FORM = compas_tno.get('form.json')
DATA_XFORM = compas_tno.get('Xform.json')
DATA_XFORCE = compas_tno.get('Xforce.json')

form = FormDiagram.from_json(DATA_FORM)
force = reciprocal_from_form(form)

animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, interval=150)
