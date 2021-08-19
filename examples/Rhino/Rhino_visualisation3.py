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
master = 'dome_polar_max_19.297'

# Put here the .json file for the optimisation
fnm = '/Users/mricardo/compas_dev/compas_tno/data/test.json'
fnm = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=5.0/pointed_crossvault_cross_fd_discr_10_min_thk_50.0.json'
fnm = '/Users/mricardo/compas_dev/me/shape_comparison/dome_polar/radial_spaced_fd/dome_polar_radial_spaced_fd_discr_[8, 20]_max_thk_19.297.json'
# fnm = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.09_discr_[4, 16]_min.json'
# fnm = '/Users/mricardo/compas_dev/me/SI_data/Amiens/pointed_crossvault/fan_fd/forms/pointed_crossvault_fan_fd_discr_14_fill_0.8_min_thk_46.0.json'
# fnm = '/Users/mricardo/compas_dev/me/minmax/dome/flower/flower_discr_8_20_min_t=50.json'
# fnm = compas_tno.get('test.json')
# fnm = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.3_discr_[4, 12]_radial_spaced_fd_straight_min.json'

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
except BaseException:
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
rs.LayerVisible(thrust_layer, False)
rs.LayerVisible(pipes_layer, True)
rs.LayerVisible(loads_layer, False)
rs.LayerVisible(pipes_color, False)
print(max(fs))
