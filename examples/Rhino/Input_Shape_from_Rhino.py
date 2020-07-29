import compas_tno
import os
import rhinoscriptsyntax as rs

from compas_rhino import select_mesh
from compas_rhino.geometry import RhinoMesh

folder = '/Users/mricardo/compas_dev/compas_tno/data/input/'

extrados_guid = select_mesh(message='Select Extrados Mesh')
extrados_RM = RhinoMesh.from_guid(extrados_guid)
extrados = extrados_RM.to_compas()
file_name = rs.GetString('Name of .json file to save', 'extrados')
extrados.to_json(os.path.join(folder, file_name + '.json'))

intrados_guid = select_mesh(message='Select Intrados Mesh')
intrados_RM = RhinoMesh.from_guid(intrados_guid)
intrados = intrados_RM.to_compas()
file_name = rs.GetString('Name of .json file to save', 'intrados')
intrados.to_json(os.path.join(folder, file_name + '.json'))

target_guid = select_mesh(message='Select Target Mesh')
if target_guid:
    target_RM = RhinoMesh.from_guid(target_guid)
    target = target_RM.to_compas()
    file_name = rs.GetString('Name of .json file to save', 'target')
    target.to_json(os.path.join(folder, file_name + '.json'))
