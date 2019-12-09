from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
# from compas_thrust.diagrams.form import overview_forces
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json
import math

# Mesh

fnm = '/Users/mricardo/compas_dev/me/minmax/pavillion/cross_fd/cross_fd_discr_20_max_t=50.json'
form = FormDiagram.from_json(fnm)
k_i = form.key_index()
i_k = form.index_key()

try:
    cracks_lb, cracks_ub = form.attributes['cracks']
except:
    cracks_lb = []
    cracks_ub = []
fs = []
radius_max = 0.30
# Usually 0.50 for Cross Form-Diagram

thrust_layer = 'Thrust' 
rs.AddLayer(thrust_layer)
rs.CurrentLayer(thrust_layer)
artist = NetworkArtist(form, layer=thrust_layer)
artist.clear_layer()
lp = 0.0
for u, v in form.edges():
    q = form.get_edge_attribute((u,v), 'q')
    l = form.edge_length(u,v)
    fs.append(q*l)
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))
    lp += q * l * l
rs.AddTextDot('{0:.1f}'.format(lp),[- 1.0, - 1.0, 0.0 ])

thrust_limits = 'Limits' 
rs.AddLayer(thrust_limits)
rs.CurrentLayer(thrust_limits)
artist = NetworkArtist(form, layer=thrust_limits)
artist.clear_layer()
lp = 0.0
tol_crack = 0.01
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    try:
        lb = form.get_vertex_attribute(key, 'lb')
        if abs(z - lb) < tol_crack:
            id = rs.AddPoint([x,y,z])
            rs.ObjectColor(id, color = (0,0,200))
    except:
        pass
    try:
        ub = form.get_vertex_attribute(key, 'ub')
        if abs(z - ub) < tol_crack:
            id = rs.AddPoint([x,y,z])
            rs.ObjectColor(id, color = (200,0,0))
    except:
        pass
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

for key in form.vertices_where({'rol_x': True}):
    try:
        rx = form.get_vertex_attribute(key, 'rx')
        sp = form.vertex_coordinates(key)
        ep = [sp[0] + rx, sp[1], sp[2]]
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(rx))
        print('Rx on key: ', key, 'Value: ', rx)
    except:
        pass

for key in form.vertices_where({'rol_y': True}):
    try:
        ry = form.get_vertex_attribute(key, 'ry')
        sp = form.vertex_coordinates(key)
        ep = [sp[0], sp[1] + ry, sp[2]]
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(ry))
        print('Ry on key: ', key, 'Value: ', ry)
    except:
        pass

pipes_layer = 'Pipes' 
rs.AddLayer(pipes_layer)
rs.CurrentLayer(pipes_layer)
artist = NetworkArtist(form, layer=pipes_layer)
artist.clear_layer()
rs.CurrentLayer(pipes_layer)
for u, v in form.edges():
    l = form.edge_length(u, v)
    q = form.get_edge_attribute((u, v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    if q > 0.001:
        id = rs.AddLine(sp, ep)
        f = math.sqrt(math.fabs(q * l))
        coef = f/max(fs)*radius_max
        pipe = rs.AddPipe(id, 0, coef)
        rs.ObjectColor(pipe, color=(255,0,0))

rs.CurrentLayer('Default')
rs.LayerVisible(target_layer, False)
rs.LayerVisible(ub_layer, False)
rs.LayerVisible(lb_layer, False)
rs.LayerVisible(thrust_layer, True)
rs.LayerVisible(reac_layer, True)
rs.LayerVisible(thrust_limits, True)
rs.LayerVisible(pipes_layer, False)