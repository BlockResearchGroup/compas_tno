
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_tna.diagrams import FormDiagram
from compas.utilities import geometric_key
from compas_rhino.artists import NetworkArtist

import rhinoscriptsyntax as rs


jsonpath = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_comp.json'
# jsonpath_complete = '/Users/mricardo/compas_dev/compas_loadpath/data/constraint/vault_comp_2.json'

# Form

Lines_txt = 'A::Lines_comp_rlx' #_complete
Symmetry_txt = 'A::Sym_comp_rlx' #_complete
Pins_txt = 'A::Pins_comp_rlx' #_complete
Dots_txt = 'A::Dots' #_complete
rollers_txt = 'Rollers'
inds_layer = 'Inds'
lb_layer = 'lb'
ub_layer = 'ub'
target_layer = '3DPoints'
dots_3D = 'Dots_3D'

guids = rs.ObjectsByLayer(Lines_txt) + rs.ObjectsByLayer(Symmetry_txt)
lines = [[rs.CurveStartPoint(i), rs.CurveEndPoint(i)] for i in guids if rs.IsCurve(i)]
form = FormDiagram.from_lines(lines, delete_boundary_face=False)

form.update_default_vertex_attributes({'is_roller': False})
form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
form.attributes['loadpath'] = 0
form.attributes['indset'] = []

gkey_key = form.gkey_key()

Loads3d = False
lb_ub = False
scale = None
rollers = False
complete = False
nsym = 2 #8
ind = False
openings = True

# Pins

for i in rs.ObjectsByLayer(Pins_txt):
    gkey = geometric_key(rs.PointCoordinates(i))
    form.set_vertex_attribute(gkey_key[gkey], 'is_fixed', True)

# Rollers

if rollers:
    for i in rs.ObjectsByLayer(rollers_txt):
        gkey = geometric_key(rs.PointCoordinates(i))
        form.set_vertex_attributes(gkey_key[gkey], {'is_roller': True})

# Loads

artist = NetworkArtist(form, layer=Dots_txt)
artist.clear_layer()
artist = NetworkArtist(form, layer=dots_3D)
artist.clear_layer()

loads = FormDiagram.from_lines(lines, delete_boundary_face=True)
pzt = 0

# Openings - if Any

if openings:
    for key in form.faces():
        if form.face_area(key) > 12.0 and form.face_area(key) < 12.2 :
            form.delete_face(key)
            print('Deleted area of face {0}'.format(key))
            break

for key in form.vertices():
    form.vertex[key]['pz'] = loads.vertex_area(key=key)
    pzt += form.vertex[key]['pz']
print('Planar load - pzt = {0}'.format(pzt))

if Loads3d == True:
   
   for i in rs.ObjectsByLayer(target_layer):
       point_target = rs.PointCoordinates(i)
       point_ground = [point_target[0],point_target[1],0.0]
       gkey = geometric_key(point_ground)
       z_target = point_target[2]
       loads.set_vertex_attribute(gkey_key[gkey], 'z', z_target)
   
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

if lb_ub == True:

    for i in rs.ObjectsByLayer(lb_layer):
       point_target = rs.PointCoordinates(i)
       point_ground = [point_target[0],point_target[1],0.0]
       gkey = geometric_key(point_ground)
       z_target = point_target[2]
       form.set_vertex_attribute(gkey_key[gkey], 'lb', z_target)
       lb_constraints += 1
       
       
    for i in rs.ObjectsByLayer(ub_layer):
       point_target = rs.PointCoordinates(i)
       point_ground = [point_target[0],point_target[1],0.0]
       gkey = geometric_key(point_ground)
       z_target = point_target[2]
       form.set_vertex_attribute(gkey_key[gkey], 'ub', z_target)
       ub_constraints += 1

    print('Got {0} lb contraints and {1} ub constraints'.format(lb_constraints,ub_constraints))

# Symmetry

for i in rs.ObjectsByLayer(Symmetry_txt):
    
    if rs.IsCurve(i):
        u = gkey_key[geometric_key(rs.CurveStartPoint(i))]
        v = gkey_key[geometric_key(rs.CurveEndPoint(i))]
        form.set_edge_attribute((u, v), name='is_symmetry', value=True)
        
    elif rs.IsPoint(i):
        u = gkey_key[geometric_key(rs.PointCoordinates(i))]
        name = rs.ObjectName(i)
        form.set_vertex_attribute(u, name='pz', value=float(name))

# Inds

if ind:
    ind = []
    for i in rs.ObjectsByLayer(inds_layer):  
        if rs.IsCurve(i):
            u = gkey_key[geometric_key(rs.CurveStartPoint(i))]
            v = gkey_key[geometric_key(rs.CurveEndPoint(i))]
            form.set_edge_attribute((u, v), name='is_ind', value=True)
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
    form.set_vertex_attribute(key, 'z', 0.0)
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

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ------------------ COMPLETE TOPOLOGY ------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


if complete:
    
    Lines_txt = 'Lines_complete' #_complete
    Symmetry_txt = 'Sym_complete' #_complete
    Pins_txt = 'Pins_complete' #_complete

    target_layer = '3DPoints_complete'

    guids = rs.ObjectsByLayer(Lines_txt) + rs.ObjectsByLayer(Symmetry_txt)
    lines = [[rs.CurveStartPoint(i), rs.CurveEndPoint(i)] for i in guids if rs.IsCurve(i)]
    form = FormDiagram.from_lines(lines, delete_boundary_face=False)

    form.update_default_vertex_attributes({'is_roller': False})
    form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form.attributes['loadpath'] = 0
    form.attributes['indset'] = []

    gkey_key = form.gkey_key()

    # Pins

    for i in rs.ObjectsByLayer(Pins_txt):
        gkey = geometric_key(rs.PointCoordinates(i))
        form.set_vertex_attribute(gkey_key[gkey], 'is_fixed', True)

    # Rollers

    if rollers:
        for i in rs.ObjectsByLayer(rollers_txt):
            gkey = geometric_key(rs.PointCoordinates(i))
            form.set_vertex_attributes(gkey_key[gkey], {'is_roller': True})

    # Loads

    loads = FormDiagram.from_lines(lines, delete_boundary_face=True)
    pzt = 0
    for key in form.vertices():
        form.vertex[key]['pz'] = loads.vertex_area(key=key)
        pzt += form.vertex[key]['pz']
    print('Planar load - pzt = {0}'.format(pzt))

    if Loads3d == True:
    
        for i in rs.ObjectsByLayer(target_layer):
            point_target = rs.PointCoordinates(i)
            point_ground = [point_target[0],point_target[1],0.0]
            gkey = geometric_key(point_ground)
            z_target = point_target[2]
            loads.set_vertex_attribute(gkey_key[gkey], 'z', z_target)
        
        pvt = 0
        for key in form.vertices():
                pz = loads.vertex_area(key=key)
                form.vertex[key]['pz'] = pz
                rs.CurrentLayer(dots_3D)
                rs.AddTextDot('{0:.2f}'.format(pz), loads.vertex_coordinates(key))
                pvt += pz

        print('3D load - pzt = {0}'.format(pvt))

    # Symmetry

    for i in rs.ObjectsByLayer(Symmetry_txt):
        
        if rs.IsCurve(i):
            u = gkey_key[geometric_key(rs.CurveStartPoint(i))]
            v = gkey_key[geometric_key(rs.CurveEndPoint(i))]
            form.set_edge_attribute((u, v), name='is_symmetry', value=True)
            
        elif rs.IsPoint(i):
            u = gkey_key[geometric_key(rs.PointCoordinates(i))]
            name = rs.ObjectName(i)
            form.set_vertex_attribute(u, name='pz', value=float(name))

    # Scale

    if scale:

        scl = scale/pvt
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
        form.set_vertex_attribute(key, 'z', 0.0)
        pzt += pz
        if pz:
            rs.AddTextDot('{0:.2f}'.format(pz), form.vertex_coordinates(key))
    print('Total load: {0}'.format(pzt))

    rs.EnableRedraw(True)
    rs.CurrentLayer('Default')

    # Save

    form.to_json(jsonpath_complete)