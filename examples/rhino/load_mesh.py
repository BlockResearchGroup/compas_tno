import compas_rhino
import compas_tno
import os
from compas_rhino.conversions import RhinoMesh
from compas_rhino.utilities import get_objects
from compas_rhino.utilities import select_mesh

while True:
    out = compas_rhino.rs.GetString('Title')
    if not out:
        break

    guid = select_mesh(message='Select the Form Diagram Mesh')

    mesh = RhinoMesh.from_guid(guid)
    compas_mesh = mesh.to_compas()

    jsonpath = '/Users/mricardo/compas_dev/me/pattern/singular/crossvault/mesh-'+ str(out) +'.json'
    #jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_square/mesh-E.json'

    compas_mesh.to_json(jsonpath)

    print('Mesh saved @', jsonpath)
    compas_rhino.rs.UnselectAllObjects()