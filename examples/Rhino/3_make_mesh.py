
from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point

import rhinoscriptsyntax as rs


guids = rs.ObjectsByLayer('Lines')
lines = [[rs.CurveStartPoint(i), rs.CurveEndPoint(i)] for i in guids if rs.IsCurve(i)]
form = FormDiagram.from_lines(lines, delete_boundary_face=True)
vertices = [form.vertex_coordinates(key) for key in form.vertices()]
faces = [form.face[key] for key in form.faces()]

rs.AddLayer('Mesh')
rs.CurrentLayer('Mesh')

mesh = rs.AddMesh(vertices,faces)