import compas_rhino
import compas_tno
import os
from compas_rhino.conversions import RhinoLine
from compas_rhino.utilities import get_objects
from compas_rhino.utilities import select_lines
from compas.datastructures import Mesh

guids = select_lines(message='Select the Lines for the Form Diagram Mesh')

lines = []
for guid in guids:
    rhinoline = RhinoLine.from_guid(guid)
    compasline = rhinoline.to_compas()
    lines.append(compasline)

compas_mesh = Mesh.from_lines(lines, delete_boundary_face=True)

#jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_square/mesh-E.json'
jsonpath = '/Users/mricardo/compas_dev/me/pattern/singular/dome/mesh-E3.json'

compas_mesh.to_json(jsonpath)

print('Mesh saved @', jsonpath)