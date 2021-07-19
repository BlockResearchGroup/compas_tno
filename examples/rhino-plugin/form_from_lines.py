from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import FormGraph
from compas_tno.rhino import Scene
import compas_rhino

scene = Scene()

guids = compas_rhino.select_lines(message='Select Form Diagram Lines')
        
lines = compas_rhino.get_line_coordinates(guids)
graph = FormGraph.from_lines(lines)

form = FormDiagram.from_lines(lines)
print(form.number_of_edges())
print(form.number_of_vertices())

scene.add(form, name="Form")
scene.save()
