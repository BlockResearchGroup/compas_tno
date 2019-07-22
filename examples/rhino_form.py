from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json


# Plot Thrust Network
i = 0

# for i in range(1,1003):

# fnm = '/Users/mricardo/compas_dev/me/discretize/02_0'+str(i)+'_complete.json'
# fnm = '/Users/mricardo/compas_dev/me/bestfit/pillow3_lp.json'
# fnm = '/Users/mricardo/compas_dev/me/bestfit/pillowRV_ind_calc.json'
# fnm = '/Users/mricardo/compas_dev/me/minmax/barrel/3D_min_ind.json'
fnm = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_calc.json'
# radical = '/Users/mricardo/compas_dev/me/loadpath/Fix/' + ST
# fnm = radical + '_lp.json'
form = FormDiagram.from_json(fnm)

lp = 0
try:
    t = form.attributes['offset']
    print('Translation of {0:.2f}'.format(t))
except:
    t = 0.0
    print('No offset!')

thrust_layer = 'A::Thrust_LP_half'
# thrust_layer = 'Thrust'
# reactions_layer = 'Reactions'

# thrust_layer = 'Thrust'
rs.AddLayer(thrust_layer)
# rs.AddLayer(reactions_layer)

# i += 1

# print(lp)

artist = NetworkArtist(form, layer=thrust_layer)
artist.clear_layer()

#form.set_vertices_attributes(form.vertices(), {'z': 0})

#for key, attr in form.vertices(True):
#    form.vertex[key]['z'] = 0.0
#    coord = form.vertex_coordinates(key)
#    point = rs.AddPoint(coord)
#    rs.ObjectName(point, str(form.vertex[key]['pz']))

rs.CurrentLayer(thrust_layer)

for uv in form.edges():
    u, v = uv
    q = form.get_edge_attribute(uv, 'q')
    print(q)
    l = form.edge_length(u,v)
    if form.get_edge_attribute((u,v), 'is_symmetry') == False:
        lp += q * l * l
        sp = form.vertex_coordinates(u)
        ep = form.vertex_coordinates(v)
        # sp[2] = form.get_vertex_attribute(u, 'target')
        # ep[2] = form.get_vertex_attribute(v, 'target')
        sp[2] += t
        ep[2] += t
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(q))
        # rs.ObjectColor(id, (255,255*(1002-i)/1002,0))

rs.AddTextDot('{0:.1f}'.format(lp),[-1.0,-1.0,0.0])

# artist = NetworkArtist(form, layer=reactions_layer)
# artist.clear_layer()
# rs.CurrentLayer(reactions_layer)

# pzt=0
# for key in form.vertices():
#     pz = form.get_vertex_attribute(key,'pz')
#     pzt+= pz
#     if form.vertex[key]['is_fixed'] is True:
#         node = form.vertex_coordinates(key)
#         node[2] += t
#         ry = form.get_vertex_attribute(key, 'ry')
#         rx = form.get_vertex_attribute(key, 'rx')
#         rs.AddTextDot('ry: {0:.1f} / rx: {0:.1f}'.format(ry,rx),node)

# print(pzt)

#rs.LayerVisible('Dots', False)
#rs.LayerVisible('Thrust', False)


# Copy Form

#rs.EnableRedraw(False)
#rs.CurrentLayer('Copy')
#rs.DeleteObjects(rs.ObjectsByLayer('Copy'))
#
#for uv in form.edges():
#    u, v = uv
#    q = form.get_edge_attribute(uv, 'q')
#    sp = form.vertex_coordinates(u)
#    ep = form.vertex_coordinates(v)
#    id = rs.AddLine(sp, ep)
#    rs.ObjectName(id, str(q))
#    
#rs.EnableRedraw(True)

# Analyse

#fnm = 'C:/compas-dev/compas_ags/data/loadpath/Ex03_nosym.json'
#form = FormDiagram.from_json(fnm)

#gkey_key = form.gkey_key()
#for guid in rs.ObjectsByLayer('Plan'):
#    q = float(rs.ObjectName(guid))
#    sp = rs.CurveStartPoint(guid)
#    ep = rs.CurveEndPoint(guid)
#    centroid = geometric_key([(sp[0] + ep[0])*0.5,(sp[1] + ep[1])*0.5,0.0])
#    print(centroid)
#    rs.AddPoint(centroid)
#    for uv in form.edges():
#        u, v = uv
#        if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) == centroid:
#            form.edge[u][v]['q'] = q
#            print('ok')
#
#form = XFunc('compas_ags.ags.loadpath3.z_from_form')(form)
#
#artist = NetworkArtist(form, layer='Analysis')
#artist.clear_layer()
#artist.draw_vertices()
#artist.draw_edges()
#
#fnm = 'C:/compas-dev/compas_ags/data/Euler/Ex03_sym_Euler_opt2-x2.json'
#form.to_json(fnm)
