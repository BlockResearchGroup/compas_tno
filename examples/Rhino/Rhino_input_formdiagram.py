from compas.datastructures import Mesh
from compas_rhino.geometry import RhinoMesh
from compas_rhino.geometry import RhinoLine
from compas_rhino.geometry import RhinoPoint
from compas.datastructures import Mesh
from compas.utilities import geometric_key
from compas.utilities import geometric_key_xy
import rhinoscriptsyntax as rs
import compas_rhino
import os

name = 'fancross'

for guid in compas_rhino.get_objects(layer='topology-'+name):
    mesh = RhinoMesh.from_guid(guid)
    compas_mesh = mesh.to_compas()
    #path = os.path.join(folder, lay + '_smooth.json')
    #compas_mesh.to_json(path)
    
# print(compas_mesh)

print(compas_mesh.number_of_edges())
from compas_tna.diagrams import FormDiagram

form = FormDiagram.from_mesh(compas_mesh)

print(form)


import math
sqrt2 = math.sqrt(2)
dist = 5.0 * sqrt2
tol = 1e-3
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    if abs(math.sqrt((x - 5.0)**2 + (y - 5.0)**2) - dist) < tol:
        form.vertex_attribute(key, 'is_fixed', True)
        print(key, x, y)
       
address = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-'+name+'/FormDiagram-'+name+'.json'

form.to_json(address)
print('saved to:', address)