import compas_tna

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tna.equilibrium import horizontal_nodal
from compas_plotters import MeshPlotter

FILE = compas_tna.get('tutorial/boundaryconditions.json')

form = FormDiagram.from_json(FILE)

from compas_tno.diagrams import FormDiagram as FormDiagramTNO
test = '/Users/mricardo/compas_dev/compas_tno/data/form_halfedge.json'
form = FormDiagramTNO.from_json(test)

force = ForceDiagram.from_formdiagram(form)

_edges = force.ordered_edges(form)
print('_edges', _edges)

# for u, v in form.edges():
#     print(form.halfedge[u][v], form.halfedge[v][u])

plotter = MeshPlotter(form, figsize=(10, 10))
plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
plotter.draw_vertices(radius=0.05)
plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_anchor': True})], radius=0.10, facecolor='000000')
plotter.draw_faces(text={key: str(key) for key in form.faces()})
plotter.show()

force.plot()

horizontal_nodal(form, force, kmax=100)

# ==============================================================================
# Visualise
# ==============================================================================

plotter = MeshPlotter(force, figsize=(12, 8), tight=True)

vertexcolor = {key: (1.0, 0.9, 0.9) for key in force.vertices() if not form.face_attribute(key, '_is_loaded')}

radius = {key: 0.05 for key in force.vertices()}
radius.update({key: 0.1 for key in force.vertices() if not form.face_attribute(key, '_is_loaded')})

plotter.draw_vertices(facecolor=vertexcolor, radius=radius)

plotter.draw_edges()

plotter.show()
