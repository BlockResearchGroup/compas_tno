from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
# from compas_thrust.diagrams.form import overview_forces
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json


# Plot Thrust Network
i = 0
# shapes = [2,4,5,7,8,10,11]
# shapes = range(2,9)


# fnm = '/Users/mricardo/compas_dev/me/bestfit/crossvault/discretize/0'+str(j)+'_0'+str(i)+'_fit_crossvault.json'
# fnm = '/Users/mricardo/compas_dev/me/bestfit/sixpartite/discretize/0'+str(j)+'_0'+str(i)+'_fit_sixpartite.json'

# Load

fnm = '/Users/mricardo/compas_dev/me/minmax/ortho/square/D_20_min_thrust+cracks.json'
form = FormDiagram.from_json(fnm)
k_i = form.key_index()

try:
    cracks_lb, cracks_ub = form.attributes['cracks']
except:
    cracks_lb = []
    cracks_ub = []

thrust_layer = 'Thrust' 
rs.AddLayer(thrust_layer)
rs.CurrentLayer(thrust_layer)
artist = NetworkArtist(form, layer=thrust_layer)
artist.clear_layer()
lp = 0.0
for u, v in form.edges():
    q = form.get_edge_attribute((u,v), 'q')
    l = form.edge_length(u,v)
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))
    lp += q * l * l
rs.AddTextDot('{0:.1f}'.format(lp),[- 1.0, - 1.0, 0.0 ])
if cracks_lb:
    for i in cracks_lb:
        rs.AddTextDot('lb',form.vertex_coordinates(k_i[i]))
if cracks_ub:
    for i in cracks_ub:
        rs.AddTextDot('ub',form.vertex_coordinates(k_i[i]))

target_layer = 'Target' 
rs.AddLayer(target_layer)
rs.CurrentLayer(target_layer)
artist = NetworkArtist(form, layer=target_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.get_edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.get_vertex_attribute(u, 'target')
    ep[2] = form.get_vertex_attribute(v, 'target')
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))


ub_layer = 'UB' 
rs.AddLayer(ub_layer)
rs.CurrentLayer(ub_layer)
artist = NetworkArtist(form, layer=ub_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.get_edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.get_vertex_attribute(u, 'ub')
    ep[2] = form.get_vertex_attribute(v, 'ub')
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))

lb_layer = 'LB' 
rs.AddLayer(lb_layer)
rs.CurrentLayer(lb_layer)
artist = NetworkArtist(form, layer=lb_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.get_edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.get_vertex_attribute(u, 'lb')
    ep[2] = form.get_vertex_attribute(v, 'lb')
    if sp[2] >= -1e-3 and ep[2] >= -1e-3:
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(q))

reac_layer = 'Reactions' 
rs.AddLayer(reac_layer)
rs.CurrentLayer(reac_layer)
artist = NetworkArtist(form, layer=reac_layer)
artist.clear_layer()
for key in form.vertices_where({'is_fixed': True}):
    node = form.vertex_coordinates(key)
    ry = form.get_vertex_attribute(key, 'ry')
    rx = form.get_vertex_attribute(key, 'rx')
    rz = form.get_vertex_attribute(key, 'rz')
    norm = (rx ** 2 + ry ** 2 + rz ** 2) ** (1/2)
    print('Rection on key: ', key)
    print(rx,ry,rz)
    if rz < 0.0 and norm > 0.0:
        sp = node
        print(sp)
        dz = rz/norm
        mult = node[2]/dz
        dz *= mult
        dx = mult* rx/norm
        dy = mult* ry/norm
        ep = [sp[0]-dx, sp[1]-dy, sp[2]-dz]
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(norm))
        rs.AddTextDot('rx: {0:.1f} / ry: {1:.1f}'.format(rx, ry), node)


rs.LayerVisible(target_layer, False)
rs.LayerVisible(ub_layer, False)
rs.LayerVisible(lb_layer, False)
rs.LayerVisible(thrust_layer, True)
rs.LayerVisible(reac_layer, True)

