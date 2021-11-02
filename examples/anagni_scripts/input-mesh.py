import compas_tno
import os
import rhinoscriptsyntax as rs

from compas_rhino import select_mesh
from compas_rhino.geometry import RhinoMesh

folder = '/Users/mricardo/compas_dev/me/anagni/'

guid = select_mesh(message='Select Mesh')
RM = RhinoMesh.from_guid(guid)
mesh = RM.to_compas()
file_name = rs.GetString('Name of .json file to save', 'top_vault_form')
address = os.path.join(folder, file_name + '.json')
mesh.to_json(address)
print('Salved to:', address)