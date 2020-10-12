from compas_tno.diagrams import FormDiagram
from compas.geometry import Point
import compas_tno

from compas_viewers.objectviewer import ObjectViewer

cracks_ub = []
cracks_lb = []

viewer = ObjectViewer()

json = compas_tno.get('test.json')
form = FormDiagram.from_json(json)

intrados = form.copy()
extrados = form.copy()

intrad = 1
extrad = 1

for key in intrados.vertices():
    lb = form.vertex_attribute(key, 'lb')
    ub = form.vertex_attribute(key, 'ub')
    x, y, z = form.vertex_coordinates(key)
    intrados.vertex_attribute(key, 'z', lb)
    extrados.vertex_attribute(key, 'z', ub)
    if abs(ub - z) < 10e-4:
        viewer.add(Point(x, y, z), name="Extrados (%s)"%extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
        extrad += 1
    elif abs(lb - z) < 10e-4:
        viewer.add(Point(x, y, z), name="Intrados (%s)"%intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
        intrad += 1

viewer.add(form, name="FormDiagram", settings={
    'color': '#FF0000',
    'edges.color': '#FF0000',
    'edges.width': 2,
    'opacity': 0.8,
    'vertices.size': 0,
    'vertices.on': False,
    'edges.on': True,
    'faces.on': False,
    })

viewer.add(intrados, name="Intrados", settings={
    'color': '#999999',
    'edges.width': 3,
    'opacity': 0.5,
    'vertices.size': 0,
    'vertices.on': False,
    'edges.on': False,
    'faces.on': True,
    })

viewer.add(extrados, name="Extrados", settings={
    'color': '#999999',
    'edges.width': 3,
    'opacity': 0.5,
    'vertices.size': 0,
    'vertices.on': False,
    'edges.on': False,
    'faces.on': True,
    })

viewer.show()
