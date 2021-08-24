import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

json_form = compas_tno.get('form.json')
json_shape = compas_tno.get('shape.json')

print(json_form)

form = FormDiagram.from_json(json_form)
shape = Shape.from_json(json_shape)


print(form)
