import compas
print(compas.__version__)
from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json


for i in range(2,8):
    j = 2
    fnm = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+str(j)+'_0'+str(i)+'_complete.json'
    fnm_base = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+str(j)+'_08_complete.json'
    
    form = FormDiagram.from_json(fnm)
    vertices = [form.vertex_coordinates(key) for key in form.vertices()]
    faces = [form.face[key] for key in form.faces()]
    mesh = rs.AddMesh(vertices,faces)
        
    # mesh_layer = str(j)+'_midsupport::Mesh_final'

    # for msh in rs.ObjectsByLayer(mesh_layer):
    #     if rs.IsMesh(msh):
    #         mesh = msh

    form.update_default_vertex_attributes({'dist': 0.0,})
    for key in form.vertices():
        if form.get_vertex_attribute(key, 'is_fixed') == False:
            pt = form.vertex_coordinates(key)
            proj3D = rs.ProjectPointToMesh(pt,mesh,(0,0,1))
            # print(proj3D)
            # print(proj3D)
            try:
                proj = [proj3D[0].X,proj3D[0].Y,proj3D[0].Z]
                dist = distance_point_point(pt, proj)
            except:
                dist = 0
                print('Warning at: {0:.3f}, {1:.3f}'.format(pt[0],pt[1]))
            # print(dist)
            form.set_vertex_attribute(key, 'dist', value=dist)

    # fnm_dist = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+str(j)+'_0'+str(i)+'_dist.json'
    # form.to_json(fnm_dist)
    # print('Form {0} saved!'.format(i))