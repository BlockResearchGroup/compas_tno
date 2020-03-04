import compas_tna

from compas_tna.diagrams import FormDiagram

FILE = compas_tna.get('tutorial/rhinomesh.obj')

print(FILE)

form = FormDiagram.from_obj(FILE)

lines = [
    [[0,0,0],[1,1,0]],
    [[0,0,0],[-5,3,0]]
]

form = FormDiagram.from_lines(lines)

corners = list(form.vertices_where({'vertex_degree': 2}))
print(corners)

for i in corners:
    form.vertex_attribute(i, 'is_anchor', True)

form.plot()

form = FormDiagram()
form.add_vertex(0,0,0)
form.add_vertex(0,0,1)
form.add_vertex(0,3,2)

