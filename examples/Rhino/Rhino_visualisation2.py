from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
from compas.utilities import i_to_rgb
# from compas_tno.diagrams.form import overview_forces
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json
import math
import compas_tno

# Give a name to your project
master = 'Projected_Straight_PX_0,30'

# Put here the .json file for the optimisation
# fnm = '/Users/mricardo/compas_dev/compas_tno/data/test.4.json'
# fnm = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.09_discr_[4, 16]_min.json'
# fnm = '/Users/mricardo/compas_dev/me/SI_data/Amiens/pointed_crossvault/fan_fd/forms/pointed_crossvault_fan_fd_discr_14_fill_0.8_min_thk_46.0.json'
# fnm = '/Users/mricardo/compas_dev/me/minmax/dome/flower/flower_discr_8_20_min_t=50.json'
# fnm = compas_tno.get('test.json')
fnm = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.3_discr_[4, 12]_radial_spaced_fd_straight_min.json'

# Parameters for the visualisation
radius_max = 0.10 # 0.175 for dome radial and 0.15 for dome flower, and 0.25 for the fan-vault
radius_circlus = 0.15 # Radius of the speres that mark when it touches upper bound and lower bound.
radius_colored_pipes = 0.10
max_length_reaction = 0.5

form = FormDiagram.from_json(fnm)
k_i = form.key_index()
i_k = form.index_key()

try:
    cracks_lb, cracks_ub = form.attributes['cracks']
except:
    cracks_lb = []
    cracks_ub = []
fs = []

thrust_layer = master + '::Thrust'
rs.AddLayer(thrust_layer)
rs.CurrentLayer(thrust_layer)
artist = NetworkArtist(form)
artist.clear_layer()
lp = 0.0
for u, v in form.edges_where({'_is_edge': True}):
    q = form.edge_attribute((u,v), 'q')
    l = form.edge_length(u,v)
    fs.append(q*l)
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))
    lp += q * l * l
rs.AddTextDot('{0:.1f}'.format(lp),[- 1.0, - 1.0, 0.0 ])

thrust_layer = master + '::Independents'
rs.AddLayer(thrust_layer)
rs.CurrentLayer(thrust_layer)
artist = NetworkArtist(form)
artist.clear_layer()
lp = 0.0
for u, v in form.edges_where({'is_ind': True}):
    q = form.edge_attribute((u,v), 'q')
    l = form.edge_length(u,v)
    fs.append(q*l)
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))
    lp += q * l * l
rs.AddTextDot('{0:.1f}'.format(lp),[- 1.0, - 1.0, 0.0 ])

thrust_limits = master + '::Limits'
rs.AddLayer(thrust_limits)
rs.CurrentLayer(thrust_limits)
artist = NetworkArtist(form, layer=thrust_limits)
artist.clear_layer()
lp = 0.0
tol_crack = 0.001
for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    try:
        lb = form.vertex_attribute(key, 'lb')
        if abs(z - lb) < tol_crack:
            id = rs.AddPoint([x,y,z])
            rs.ObjectColor(id, color = (0,0,255))
            id = rs.AddSphere([x,y,z], radius_circlus)
            rs.ObjectColor(id, color = (0,0,255))
    except:
        pass
    try:
        ub = form.vertex_attribute(key, 'ub')
        if abs(z - ub) < tol_crack:
            id = rs.AddPoint([x,y,z])
            rs.ObjectColor(id, color = (0,128,0))
            id = rs.AddSphere([x,y,z], radius_circlus)
            rs.ObjectColor(id, color = (0,128,0))
    except:
        pass
if cracks_lb:
    for i in cracks_lb:
        rs.AddTextDot('lb',form.vertex_coordinates(k_i[i]))
if cracks_ub:
    for i in cracks_ub:
        rs.AddTextDot('ub',form.vertex_coordinates(k_i[i]))

target_layer = master + '::Target'
rs.AddLayer(target_layer)
rs.CurrentLayer(target_layer)
artist = NetworkArtist(form, layer=target_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.vertex_attribute(u, 'target')
    ep[2] = form.vertex_attribute(v, 'target')
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))

loads_layer = master + '::Loads'
rs.AddLayer(loads_layer)
rs.CurrentLayer(loads_layer)
artist = NetworkArtist(form, layer=loads_layer)
artist.clear_layer()
for key in form.vertices():
    node = form.vertex_coordinates(key)
    pz = form.vertex_attribute(key, 'pz')
    rs.AddTextDot('({0:.1f})'.format(pz), node, )

ub_layer = master + '::UB'
rs.AddLayer(ub_layer)
rs.CurrentLayer(ub_layer)
artist = NetworkArtist(form, layer=ub_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.vertex_attribute(u, 'ub')
    ep[2] = form.vertex_attribute(v, 'ub')
    id = rs.AddLine(sp, ep)
    rs.ObjectName(id, str(q))

lb_layer = master + '::LB'
rs.AddLayer(lb_layer)
rs.CurrentLayer(lb_layer)
artist = NetworkArtist(form, layer=lb_layer)
artist.clear_layer()
for u, v in form.edges():
    q = form.edge_attribute((u,v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    sp[2] = form.vertex_attribute(u, 'lb')
    ep[2] = form.vertex_attribute(v, 'lb')
    if sp[2] >= -1e-3 and ep[2] >= -1e-3:
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(q))

reac_layer = master + '::Reactions'
reac_val = master + '::Reactions-val'
rs.AddLayer(reac_layer)
rs.AddLayer(reac_val)
rs.CurrentLayer(reac_layer)
artist = NetworkArtist(form, layer=reac_layer)
artist.clear_layer()
for key in form.vertices_where({'is_fixed': True}):
    node = form.vertex_coordinates(key)
    ry = form.vertex_attribute(key, '_ry')
    rx = form.vertex_attribute(key, '_rx')
    rz = form.vertex_attribute(key, '_rz')
    norm = (rx ** 2 + ry ** 2 + rz ** 2) ** (1/2)
    if rz < 0.0 and norm > 0.0:
        sp = node
        dz = rz/norm
        mult = node[2]/dz
        dz *= mult
        dx = mult* rx/norm
        dy = mult* ry/norm
        ep = [sp[0]-dx, sp[1]-dy, sp[2]-dz]
        #id = rs.AddLine(sp, ep)
        #rs.ObjectName(id, str(norm))
        rs.CurrentLayer(reac_val)
        rs.AddTextDot('({0:.1f};{1:.1f})'.format(rx, ry), node, )
        rs.CurrentLayer(reac_layer)

for key in form.vertices_where({'rol_x': True}):
    try:
        rx = form.vertex_attribute(key, '_rx')
        sp = form.vertex_coordinates(key)
        ep = [sp[0] + rx, sp[1], sp[2]]
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(rx))
        print('Rx on key: ', key, 'Value: ', rx)
    except:
        pass

for key in form.vertices_where({'rol_y': True}):
    try:
        ry = form.vertex_attribute(key, '_ry')
        sp = form.vertex_coordinates(key)
        ep = [sp[0], sp[1] + ry, sp[2]]
        id = rs.AddLine(sp, ep)
        rs.ObjectName(id, str(ry))
        print('Ry on key: ', key, 'Value: ', ry)
    except:
        pass

rs.LayerVisible(reac_val, False)
pipes_layer = master + '::Pipes'
rs.AddLayer(pipes_layer)
rs.CurrentLayer(pipes_layer)
artist = NetworkArtist(form, layer=pipes_layer)
artist.clear_layer()
rs.CurrentLayer(pipes_layer)
for u, v in form.edges_where({'_is_edge': True}):
    l = form.edge_length(u, v)
    q = form.edge_attribute((u, v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    if q > 0.001:
        id = rs.AddLine(sp, ep)
        f = math.fabs(q * l)
        coef = f/max(fs)*radius_max
        pipe = rs.AddPipe(id, 0, coef)
        rs.ObjectColor(pipe, color=(255,0,0))
        rs.DeleteObject(id)
pipes_color = master + '::Pipes-Color'
rs.AddLayer(pipes_color)
rs.CurrentLayer(pipes_color)
artist = NetworkArtist(form, layer=pipes_color)
artist.clear_layer()
rs.CurrentLayer(pipes_color)
for u, v in form.edges_where({'_is_edge': True}):
    l = form.edge_length(u, v)
    q = form.edge_attribute((u, v), 'q')
    sp = form.vertex_coordinates(u)
    ep = form.vertex_coordinates(v)
    if q > 0.001:
        id = rs.AddLine(sp, ep)
        f = math.fabs(q * l)
        coef = f/max(fs)
        r, g, b = i_to_rgb(coef)
        rs.ObjectColor(id, color=(r,g,b))
        pipe = rs.AddPipe(id, 0, radius_colored_pipes)
        rs.ObjectColor(pipe, color=(r,g,b))
        rs.DeleteObject(id)
rs.CurrentLayer('Default')
rs.LayerVisible(target_layer, False)
rs.LayerVisible(ub_layer, False)
rs.LayerVisible(lb_layer, False)
rs.LayerVisible(thrust_layer, False)
rs.LayerVisible(reac_layer, False)
rs.LayerVisible(thrust_limits, True)
rs.LayerVisible(pipes_layer, True)
rs.LayerVisible(loads_layer, False)
rs.LayerVisible(pipes_color, False)
print(max(fs))
