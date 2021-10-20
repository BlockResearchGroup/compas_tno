import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer

path = compas_tno.get('')
address = os.path.join(path, 'form.json')

form = FormDiagram.from_json(address)

for vertex in form.vertices():
    x, y, z = form.vertex_coordinates(vertex)
    ub = form.vertex_attribute(vertex, 'ub')
    lb = form.vertex_attribute(vertex, 'lb')
    tub = form.vertex_attribute(vertex, 'tub')
    tlb = form.vertex_attribute(vertex, 'tlb')
    print('Vertex #:', vertex, ' Coordinates: ', x, y, z, 'Bounds:', ub, lb, 'ts:', tub, tlb)

for edge in form.edges():
    u, v = edge
    q = form.edge_attribute(edge, 'q')
    length = form.edge_length(u, v)
    f = q * length
    print('Edge #:', edge, ' force: ', f)


view = Viewer(form)
view.show_solution()
