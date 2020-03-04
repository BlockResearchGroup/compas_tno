
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_tna.diagrams import FormDiagram
from compas.utilities import geometric_key
# from compas_rhino.artists import NetworkArtist

import rhinoscriptsyntax as rs

# jsonpath = '/Users/mricardo/compas_dev/me/minmax/2D_arch/01_joints.json'
jsonpath = '/Users/mricardo/compas_dev/me/minmax/dome/par-diag/par-diag_discr_8_16.json'



# for i in range[(2,9)]:
j = 2
i = 6
# jsonpath_complete = '/Users/mricardo/compas_dev/compas_loadpath/data/constraint/vault_comp_2.json'
# jsonpath = '/Users/mricardo/compas_dev/me/loadpath/Fix/nosym/0'+str(j)+'_0'+str(i)+'_t_60.json'
# jsonpath = '/Users/mricardo/compas_dev/me/loadpath/freeform/test.json'
# jsonpath = '/Users/mricardo/compas_dev/me/convex/4bars/diagram.json'


# Form

Lines_txt = 'Lines'#_0' + str(j) #_complete '0' + str(j) + '_0' + str(i)
Lines_diag = 'Lines-diag'#_0' + str(j) #_complete '0' + str(j) + '_0' + str(i)
Symmetry_txt = 'Sym' #_complete
Pins_txt = 'Pins' #_complete
Dots_txt = 'Dots' #_complete
rollers_txt = 'Rollers'
inds_layer = 'Inds'
lb_layer = 'lb'
ub_layer = 'ub'
target_layer = 'target'
dots_3D = 'Dots_3D'
buttress_Layer = 'Buttress'
joints_layer = 'Joint_Segments'

guids_lines = rs.ObjectsByLayer(Lines_txt)
guids = rs.ObjectsByLayer(Lines_txt) + rs.ObjectsByLayer(Lines_diag)
lines = [[rs.CurveStartPoint(i), rs.CurveEndPoint(i)] for i in guids if rs.IsCurve(i)]
lines_main = [[rs.CurveStartPoint(i), rs.CurveEndPoint(i)] for i in guids_lines if rs.IsCurve(i)]
form = FormDiagram.from_lines(lines, delete_boundary_face=True) # AAAALWAYS CHECK IT

form.update_default_vertex_attributes({'is_roller': False})
form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
form.attributes['loadpath'] = 0
form.attributes['indset'] = []

gkey_key = form.gkey_key()

Loads3d = False # Change
lb_ub = False # Change
target = False # Change
scale = False
rollers = False
writepz = False
nsym = 1 #8
ind = False
openings = False
buttress = False
joints = False

# Pins

for i in rs.ObjectsByLayer(Pins_txt):
    gkey = geometric_key(rs.PointCoordinates(i))
    form.vertex_attribute(gkey_key[gkey], 'is_fixed', True)

# Rollers

if rollers:
    for i in rs.ObjectsByLayer(rollers_txt):
        gkey = geometric_key(rs.PointCoordinates(i))
        form.vertex_attributes(gkey_key[gkey], {'is_roller': True})

# Loads

# artist = NetworkArtist(form, layer=Dots_txt)
# artist.clear_layer()
# artist = NetworkArtist(form, layer=dots_3D)
# artist.clear_layer()

loads = FormDiagram.from_lines(lines_main, delete_boundary_face=True)
pzt = 0

# Openings - if Any

if openings:
    for key in form.faces():
        if form.face_area(key) > openings - 1.0 and form.face_area(key) < openings + 1.0:
            form.delete_face(key)
            print('Deleted area of face {0}'.format(key))
            break
    for key in loads.faces():
        if loads.face_area(key) > openings - 1.0 and loads.face_area(key) < openings + 1.0:
            loads.delete_face(key)
            print('Deleted area of face {0}'.format(key))
            break


for key in form.vertices():
    form.vertex[key]['pz'] = 1.0# loads.vertex_area(key=key)
    pzt += form.vertex[key]['pz']
print('Planar load - pzt = {0}'.format(pzt))

# Butress - if Any

if buttress:
    for i in rs.ObjectsByLayer(buttress_Layer):
        sp = rs.CurveStartPoint(i)
        ep = rs.CurveEndPoint(i)
        if rs.IsCurve(i):
            try:
                key = gkey_key[geometric_key(sp)]
                b = ep - sp
            except:
                key = gkey_key[geometric_key(ep)]
                b = sp - ep
            form.vertex_attribute(key, name = 'b', value = [b[0], b[1]])
            print('Butress for key: {0}'.format(key))
            print([b[0], b[1]])
            print(form.vertex_coordinates(key))

# Joints if desired

if joints:
    from compas.geometry import is_point_on_segment
    joints_list = []
    for i in rs.ObjectsByLayer(joints_layer):
        sp = rs.CurveStartPoint(i)
        ep = rs.CurveEndPoint(i)
        mid = rs.CurveMidPoint(i)
        mid.Z = 0.0
        # mid = 0.5 * [float(sp.X) + float(ep.X), float(sp.Y) + float(ep.Y), 0.0]
        for u,v in form.edges():
            edge = [form.vertex_coordinates(u),form.vertex_coordinates(v)]
            if is_point_on_segment(mid,edge,tol=1e-6) == True:
                joints_list.append([[sp.X,sp.Y,sp.Z],[ep.X,ep.Y,ep.Z],(u,v)])
    form.attributes['joints'] = joints_list
    print(joints_list)


if Loads3d == True:

    for i in rs.ObjectsByLayer(target_layer):
        point_target = rs.PointCoordinates(i)
        point_ground = [point_target[0],point_target[1],0.0]
        gkey = geometric_key(point_ground)
        z_target = point_target[2]
        loads.vertex_attribute(gkey_key[gkey], 'z', z_target)

    pvt = 0
    for key in form.vertices():
            pz = loads.vertex_area(key=key)
            form.vertex[key]['pz'] = pz
            rs.CurrentLayer(dots_3D)
            rs.AddTextDot('{0:.2f}'.format(pz), loads.vertex_coordinates(key))
            pvt += pz

    print('3D load - pzt = {0}'.format(pvt))

# Constraints

lb_constraints = 0
ub_constraints = 0
targets = 0

if lb_ub == True:

    for i in rs.ObjectsByLayer(lb_layer):
        point_target = rs.PointCoordinates(i)
        point_ground = [point_target[0],point_target[1],0.0]
        gkey = geometric_key(point_ground)
        z_target = point_target[2]
        form.vertex_attribute(gkey_key[gkey], 'lb', z_target)
        lb_constraints += 1


    for i in rs.ObjectsByLayer(ub_layer):
        point_target = rs.PointCoordinates(i)
        point_ground = [point_target[0],point_target[1],0.0]
        gkey = geometric_key(point_ground)
        z_target = point_target[2]
        form.vertex_attribute(gkey_key[gkey], 'ub', z_target)
        ub_constraints += 1

if target == True:

    for i in rs.ObjectsByLayer(target_layer):
        point_target = rs.PointCoordinates(i)
        point_ground = [point_target[0],point_target[1],0.0]
        gkey = geometric_key(point_ground)
        z_target = point_target[2]
        form.vertex_attribute(gkey_key[gkey], 'target', z_target)
        targets += 1

    print('Got {0} lb contraints {1} ub and {2} target constraints'.format(lb_constraints,ub_constraints,targets))

# Symmetry

for i in rs.ObjectsByLayer(Symmetry_txt):

    if rs.IsCurve(i):
        u = gkey_key[geometric_key(rs.CurveStartPoint(i))]
        v = gkey_key[geometric_key(rs.CurveEndPoint(i))]
        form.edge_attribute((u, v), name='is_symmetry', value=True)
    if writepz == False:
        if rs.IsPoint(i):
            u = gkey_key[geometric_key(rs.PointCoordinates(i))]
            name = rs.ObjectName(i)
            form.vertex_attribute(u, name='pz', value=float(name))
    else:
        if rs.IsPoint(i):
            u = gkey_key[geometric_key(rs.PointCoordinates(i))]
            pz = form.vertex_attribute(u,'pz')
            rs.ObjectName(i, str(pz))

# Inds

if ind:
    ind = []
    for i in rs.ObjectsByLayer(inds_layer):
        if rs.IsCurve(i):
            u = gkey_key[geometric_key(rs.CurveStartPoint(i))]
            v = gkey_key[geometric_key(rs.CurveEndPoint(i))]
            form.edge_attribute((u, v), name='is_ind', value=True)
            key = form.edge_midpoint(u, v)[:2] + [0.0]
            print(key)
            ind.append(geometric_key(key))
    form.attributes['indset'] = ind

# Scale

if scale:
    scl = scale/nsym/pvt
    for key in form.vertices():
        form.vertex[key]['pz'] *= scl

    print('Loads Scaled to: {0} / Scale Factor: {1}'.format(scale/nsym,scl))

# TextDots

rs.EnableRedraw(False)
rs.DeleteObjects(rs.ObjectsByLayer(Dots_txt))
rs.CurrentLayer(Dots_txt)

pzt = 0
for key in form.vertices():
    pz = form.vertex[key].get('pz', 0)
    form.vertex_attribute(key, 'z', 0.0)
    pzt += pz
    if pz:
        rs.AddTextDot('{0:.2f}'.format(pz), form.vertex_coordinates(key))
print('Total load: {0}'.format(pzt))

rs.EnableRedraw(True)
rs.CurrentLayer('Default')

# Save

print(form.number_of_edges())
form.to_json(jsonpath)
print(jsonpath)
