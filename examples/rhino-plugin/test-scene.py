from compas_tno.rhino import FormObject
from compas_tno.rhino import Scene

from compas_tno.diagrams import FormDiagram

scene = Scene()

JSON = '/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=0/crossvault_fan_fd_discr_14_deg=0_min_thk_50.0.json'
form = FormDiagram.from_json(JSON)

guid = scene.add(form, name="Form", layer="TNO::FormDiagram")

print(guid)