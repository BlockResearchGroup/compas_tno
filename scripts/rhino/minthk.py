from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram

form_ad = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk.json'
force_ad = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_force.json'

form = FormDiagram.from_json(form_ad)
force = ForceDiagram.from_json(force_ad)

