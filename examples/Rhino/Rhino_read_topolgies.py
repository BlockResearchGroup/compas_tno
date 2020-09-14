from compas.datastructures import Mesh
from compas_rhino.geometry import RhinoMesh
import compas_rhino
import os

folder = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-domical'
layers = ['Mesh']

for lay in layers:
    layer_dense = lay+'::mesh'
    for guid in compas_rhino.get_objects(layer=layer_dense):
        mesh = RhinoMesh.from_geometry(guid)
        compas_mesh = mesh.to_compas()
        path = os.path.join(folder, lay + '-domical.json')
        compas_mesh.to_json(path)

#    for guid in compas_rhino.get_objects(layer=layer_smooth):
#        mesh = RhinoMesh.from_guid(guid)
#        compas_mesh = mesh.to_compas()
#        path = os.path.join(folder, lay + '_smooth.json')
#        compas_mesh.to_json(path)
#        
print('end')