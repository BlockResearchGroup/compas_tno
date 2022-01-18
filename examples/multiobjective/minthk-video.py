from compas_tna.diagrams.formdiagram import FormDiagram
import compas_tno
from compas_tno.viewers import animation_from_optimisation

# ---- MAKE THE VIDEO ----

DATA_FORM = compas_tno.get('form.json')
DATA_XFORM = compas_tno.get('Xform.json')

form = FormDiagram.from_json(DATA_FORM)
animation_from_optimisation(form, DATA_XFORM, interval=150)
