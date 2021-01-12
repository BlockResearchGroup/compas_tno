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

folder = '/Users/mricardo/compas_dev/me/freeform/gene/'
layers = ['RV2::Lines']
# layers = ['RV2::ThrustDiagram']
lines = []

for lay in layers:
    for guid in compas_rhino.get_objects(layer=lay):
        line = RhinoLine.from_guid(guid)
        end = rs.CurveEndPoint(guid)
        start = rs.CurveStartPoint(guid)
        l = [[start.X, start.Y, start.Z],[end.X,end.Y,end.Z]]
        lines.append(l)

mesh = Mesh.from_lines(lines, delete_boundary_face=True)
print(mesh.number_of_edges())
print(mesh.number_of_faces())
print(mesh.number_of_vertices())
gkey_key = mesh.gkey_key()

layers = ['RV2::Supports']

supports = []
i=0
for lay in layers:
    for guid in compas_rhino.get_objects(layer=lay):
        i+=1
        pt = RhinoPoint.from_guid(guid)
        point = [pt.x, pt.y, pt.z]
        key = gkey_key[geometric_key(point)]
        mesh.vertex_attribute(key, 'is_fixed', True)
        mesh.vertex_attribute(key, 'target', 0.0)


print('length sup', i)

layers = ['RV2::Target']

gkey_key_xy = {geometric_key_xy(mesh.vertex_coordinates(key)): key for key in mesh.vertices()}

target = []
i=0
for lay in layers:
    for guid in compas_rhino.get_objects(layer=lay):
        i+=1
        pt = RhinoPoint.from_guid(guid)
        point = [pt.x, pt.y, pt.z]
        key = gkey_key_xy[geometric_key_xy(point)]
        mesh.vertex_attribute(key, 'target', pt.z)
        if key == 248:
            print(point)

print('length target', i)

path = os.path.join(folder, 'mesh.json')
mesh.to_json(path)    

#    for guid in compas_rhino.get_objects(layer=layer_smooth):
#        mesh = RhinoMesh.from_guid(guid)
#        compas_mesh = mesh.to_compas()
#        path = os.path.join(folder, lay + '_smooth.json')
#        compas_mesh.to_json(path)
#        
print('end')