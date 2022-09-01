import compas_rhino
import compas_tno
import os
import compas_rhino
from compas_rhino.conversions import RhinoLine
from compas_rhino.utilities import get_objects
from compas_rhino.utilities import select_lines
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.geometry import distance_point_point_xy

guids = select_lines(message='Select the Lines for the Form Diagram Mesh')

lines = []
for guid in guids:
    rhinoline = RhinoLine.from_guid(guid)
    compasline = rhinoline.to_compas()
    lines.append(compasline)
    
compas_rhino.rs.UnselectObjects(guids)

guids = select_lines(message='Select the Independent edges in the Form Diagram')

independent_lines = []
for guid in guids:
    rhinoline = RhinoLine.from_guid(guid)
    compasline = rhinoline.to_compas()
    independent_lines.append(compasline)

compas_rhino.rs.UnselectObjects(guids)

lines.sort()

compas_mesh = Mesh.from_lines(lines, delete_boundary_face=True)

ni = 4
ns = 2

title = 'circular_sigularity_ni_{}_ns_{}.json'.format(ni, ns)

# jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_square/mesh-E.json'
# jsonpath = '/Users/mricardo/compas_dev/me/anagni/meshes/CISM/mesh-B7-maxload.json'
# jsonpath = '/Users/mricardo/compas_dev/compas_tno/data/form-hor_load.json'

jsonpath= '/Users/mricardo/compas_dev/me/pattern/dome/split_support/parametric/'+title+'.json'
# jsonpath= '/Users/mricardo/compas_dev/me/pattern/singular/dome/mesh-D3-diag.json'

# jsonpath= '/Users/mricardo/compas_dev/me/pattern/parametric/form_lambd_0.5_from_rhino.json'

compas_mesh.to_json(jsonpath)

print(compas_mesh)

print('Mesh saved @', jsonpath)

from compas_tno.diagrams import FormDiagram

form = FormDiagram.from_mesh(compas_mesh)

form.set_boundary_supports()
form.delete_boundary_edges()

indset = [line.midpoint for line in independent_lines]
tol_old_ind = 1e-3
inds = []
for u, v in form.edges_where({'_is_edge': True}):
    edgemid = Point(*(form.edge_midpoint(u, v)[:2] + [0]))
    for pt in indset:
        if distance_point_point_xy(edgemid, pt) < tol_old_ind:
            inds.append((u, v))

form.assign_inds(inds)

formpath= '/Users/mricardo/compas_dev/me/pattern/dome/split_support/parametric/form-'+title+'.json'

form.to_json(formpath)

print('Form saved @', formpath)
