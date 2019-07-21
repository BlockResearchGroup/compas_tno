

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from compas.utilities import geometric_key

from compas.geometry.transformations.transformations import mirror_point_line
from compas.geometry.transformations.transformations import rotate_points
from compas.geometry import closest_point_on_line
from compas.geometry import midpoint_line_xy
from compas.geometry import matrix_from_axis_and_angle

from compas.geometry.distance import distance_point_point_xy
from numpy import argmin
from compas_thrust.algorithms.equilibrium import z_from_form

from compas_plotters import MeshPlotter

import math

from compas_thrust.plotters.plotters import plot_form

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'replicate',
    'overwiew_forces',
    'check_constraints',
    'replicate_constraints',
    'interp_surf',
    'null_edges',
    'create_sym',
]

def replicate(form,file, plot=None):

    """ Copy 1/8 form-diagram in a complete structure.

    Parameters
    ----------
    form : obj
        FormDiagram 1/8.

    form : obj
        FormDiagram complete.

    Returns
    -------
    None

    """
    from compas.utilities import reverse_geometric_key
    form_ = FormDiagram.from_json(file)

    if plot:
        plot_form(form_, radius=0.05).show()

    q_i = {}
    mid_ind = []

    for u, v in form.edges():
        if form.get_edge_attribute((u,v),'is_symmetry') == False:
            sp, ep = form.vertex_coordinates(u), form.vertex_coordinates(v)
            line_edge = [sp, ep]
            mp = midpoint_line_xy(line_edge)

            if mp[0] < 0.001:
                mp[0] = 0.0
            if mp[1] < 0.001:
                mp[1] = 0.0

            gkey_mid = geometric_key(mp)
            q_i[gkey_mid] = form.get_edge_attribute((u,v),'q')
            line_sym = ([0.0,10.0,0.0],[10.0,0.0,0.0])

            mirror_diag = mirror_point_line(mp,line_sym)

            if mirror_diag[0] < 0.001:
                mirror_diag[0] = 0.0
            if mirror_diag[1] < 0.001:
                mirror_diag[1] = 0.0


            gkey_mirror = geometric_key(mirror_diag)

            if gkey_mirror == gkey_mid:
                q_i[gkey_mid] = 2 * q_i[gkey_mid]
            elif round(mp[0], 3) == 5.000:
                q_i[gkey_mid] = 2 * q_i[gkey_mid]
                q_i[gkey_mirror] = q_i[gkey_mid]
            else:
                q_i[gkey_mirror] = q_i[gkey_mid]

            for j in range(2, 7, 2):
                [rot,diag_rot] = rotate_points([mp,mirror_diag], j/2.0*math.pi/2, axis=[0.0,0.0,1.0],origin=[5.0,5.0,0])

                if rot[0] < 0.001:
                    rot[0] = 0.0
                if rot[1] < 0.001:
                    rot[1] = 0.0

                if diag_rot[0] < 0.001:
                    diag_rot[0] = 0.0
                if diag_rot[1] < 0.001:
                    diag_rot[1] = 0.0

                gkey_rot = geometric_key(rot)
                gkey_diag_rot = geometric_key(diag_rot)

                q_i[gkey_rot] = q_i[gkey_mid]
                q_i[gkey_diag_rot] = q_i[gkey_mirror]

        if form.get_edge_attribute((u,v),'is_ind') == True:
            sp, ep = form.vertex_coordinates(u), form.vertex_coordinates(v)
            line_edge = [sp, ep]
            mp = midpoint_line_xy(line_edge)

            if mp[0] < 0.001:
                mp[0] = 0.0
            if mp[1] < 0.001:
                mp[1] = 0.0

            gkey_mid = geometric_key(mp)
            mid_ind.append(gkey_mid)

    real_points = []
    for key in q_i:
        real_points.append(reverse_geometric_key(key))

    for u, v in form_.edges():
        sp, ep = form_.vertex_coordinates(u), form_.vertex_coordinates(v)
        line_edge = [sp, ep]
        mp = midpoint_line_xy(line_edge)

        if mp[0] < 0.001:
                mp[0] = 0.0
        if mp[1] < 0.001:
            mp[1] = 0.0

        gkey_mid = geometric_key(mp)
        try:
            form_.set_edge_attribute((u, v), name='q', value = q_i[gkey_mid])
        except:
            dist = []
            for point in real_points:
                dist.append(distance_point_point_xy(point,mp))
            point_i = argmin(dist)
            print(point_i)
            gkey_appx = geometric_key(real_points[point_i])
            form_.set_edge_attribute((u, v), name='q', value = q_i[gkey_appx])

        if gkey_mid in mid_ind:
            form_.set_edge_attribute((u, v), name='is_ind', value = True)
        else:
            form_.set_edge_attribute((u, v), name='is_ind', value = False)


    # plot_form(form_, radius=0.05).show()

    try:
        t = form.attributes['offset']
    except:
        t = 0.0

    form_.attributes['offset'] = t

    form_ = z_from_form(form_)

    ql2 = 0.0
    for u, v in form_.edges():
        ql2 += form_.get_edge_attribute((u,v),'q') * form_.edge_length(u,v) ** 2

    form_.attributes['loadpath'] = ql2
    oveview_forces(form_)

    if plot:
        plot_form(form_, radius=0.05).show()

    return form_

def oveview_forces(form):

    f = []
    q = []
    z = []
    pz = 0

    lp=0

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True and form.get_edge_attribute((u,v),'is_symmetry') is False:
            qi = form.get_edge_attribute((u,v),'q')
            li = form.edge_length(u,v)
            lp += qi*li**2
            q.append(qi)
            f.append(qi*li)

    print('='*20)
    print('Overview on forces:')

    print('q: {0:.3f} : {1:.3f}'.format(float(min(q)), float(max(q))))
    print('f: {0:.3f} : {1:.3f}'.format(float(min(f)), float(max(f))))
    for key in form.vertices():
        z.append(form.get_vertex_attribute(key,'z'))
        pz += form.get_vertex_attribute(key,'pz')
    print('z: {0:.3f} : {1:.3f}'.format(float(min(z)), float(max(z))))
    print('pz: {0:.3f}'.format(pz))
    print('lp: {0:.3f}'.format(lp))

    return

def check_constraints(form, show=False):

    try:
        t = form.attributes['offset']
    except:
        t = 0.0
    outside = {}
    penalty = 0

    for key, vertex in form.vertex.items():
        z = form.vertex_coordinates(key)[2] + t
        if vertex.get('lb', None):
            lb = vertex['lb']
            if z < lb:
                outside[key] = lb - z
                penalty += (abs(outside[key])+4)**(4)
        if vertex.get('ub', None):
            ub = vertex['ub']
            if z > ub:
                outside[key] = z - ub
                penalty += (abs(outside[key])+4)**(4)

    print('The penalty in the constraints is {0:.3f}'.format(penalty))

    if show:
        plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
        plotter.draw_vertices(text=outside)
        plotter.draw_edges()
        plotter.show()

    return penalty

def replicate_contraints(file, file_constraint):

    form = FormDiagram.from_json(file)
    form_ = FormDiagram.from_json(file_constraint)
    gkey_planar = {}

    for key_real in form.vertices():
        coord = form.vertex_coordinates(key_real)
        gkey_proj = geometric_key([coord[0],coord[1],0.0])
        gkey_planar[gkey_proj] = key_real

    for key in form_.vertices():
        target = form_.vertex[key].get('target', 0.0)
        if target < 10**(-4):
            target = 0.00
        gkey = geometric_key([form_.vertex_coordinates(key)[0],form_.vertex_coordinates(key)[1], 0.0])
        form.set_vertex_attribute(gkey_planar[gkey], 'target', target)

    return form

def interp_surf(form):

    x = []
    y = []
    s = []

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            x.append(vertex.get('x'))
            y.append(vertex.get('y'))
            s.append(vertex.get('target'))

    from scipy import interpolate

    surf = interpolate.interp2d(x, y, s, kind = 'linear')

    return surf

def null_edges(form, plot=False):

    null_edges = []
    all_edges = []

    for u, v in form.edges():
        if form.get_edge_attribute((u,v), 'is_external') == False and form.get_edge_attribute((u,v), 'is_edge') == True:
            activ = 0
            coord_u = form.vertex_coordinates(u)
            coord_v = form.vertex_coordinates(v)
            ux = round(coord_u[0],3)
            uy = round(coord_u[1],3)
            vx = round(coord_v[0],3)
            vy = round(coord_v[1],3)
            mid_x, mid_y, _ = form.edge_midpoint(u,v)
            if uy == vy and ((uy is not 10.0 and vy is not 10.0) or (uy is not 0.0 and vy is not 0.0)):
                if (mid_y > mid_x and mid_y < 10 - mid_x) or (mid_y < mid_x and mid_y > 10 - mid_x):
                    if uy == 5.0 and vy == 5.0 and ux > 0.01 and vx > 0.01 and ux < 9.99 and vx < 9.99: # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u,v))
                        activ += 1
            if ux == vx and ((ux is not 10.0 and vx is not 10.0) or (ux is not 0.0 and vx is not 0.0)):
                if (mid_y > mid_x and mid_y > 10 - mid_x) or (mid_y < mid_x and mid_y < 10 - mid_x):
                    if ux == 5.0 and vx == 5.0 and uy > 0.01 and vy > 0.01 and uy < 9.99 and vy < 9.99: # Special for TOP 2
                        pass
                    else:
                        null_edges.append((u,v))
                        activ += 1
            if activ == 0:
                all_edges.append((u,v))

    if plot:
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges(all_edges)
        plotter.draw_edges(null_edges, color='ff0000')
        plotter.show()

    return null_edges

def create_sym(form):

    plot_form(form).show()

    lines = []
    symmetry = []
    pins = []
    loads = {}
    tol = 0.001

    for u,v in form.edges():
        u_coord = form.vertex_coordinates(u)
        v_coord = form.vertex_coordinates(v)
        if u_coord[1] >= 10 - u_coord[0] - tol and u_coord[0] <= 5.0 + tol and v_coord[1] >= 10 - v_coord[0] - tol and v_coord[0] <= 5.0 + tol :
            lines.append([u_coord,v_coord])

    for key in form.vertices():
        coord = form.vertex_coordinates(key)
        gkey = geometric_key(coord)
        loads[gkey] = form.get_vertex_attribute(key, 'pz')

        if 5.0 - tol <= coord[0] <= 5.0 + tol and coord[1] > 5.0 - tol: # Line Vertical at x = 5.0 and y > 5.0
            loads[gkey] = 0.5 * loads[gkey]
            if coord[1] < 10.0 - tol:  # Only if continuously supported
                ext_pt = [coord[0]+1.0, coord[1], coord[2]]
                symmetry.append([coord,ext_pt])
                pins.append(geometric_key(ext_pt))
                loads[geometric_key(ext_pt)] = 0.0

        if 10 - coord[0] - tol <= coord[1] <= 10 - coord[0] + tol and 10.0 - tol > coord[1] > 5.0 + tol: # Line Diagonal at y = x - 10 without central and corner
            loads[gkey] = 0.5 * loads[gkey]
            ext_pt = [coord[0]-sqrt(2)/2,coord[1]-sqrt(2)/2,coord[2]]
            symmetry.append([coord,ext_pt])
            pins.append(geometric_key(ext_pt))
            loads[geometric_key(ext_pt)] = 0.0

        if 5.0 - tol <= coord[0] <= 5.0 + tol and 5.0 - tol <= coord[1] <= 5.0 + tol: # Central Point
            loads[gkey] = 0.25 * loads[gkey]
            ext_pt = [coord[0],coord[1]-1.0,coord[2]]
            symmetry.append([coord,ext_pt])
            pins.append(geometric_key(ext_pt))
            loads[geometric_key(ext_pt)] = 0.0

        if form.get_vertex_attribute(key, 'is_fixed') == True and coord[1] > 10.0 - tol and tol <= coord[0] <= 5.0 + tol: # Superior line without corner point
            pins.append(geometric_key(coord))

        if form.get_vertex_attribute(key, 'is_fixed') == True and coord[1] > 10.0 - tol and -tol <= coord[0] <= tol: # Corner Point
            pins.append(geometric_key(coord))
            loads[gkey] = 0.5 * loads[gkey]

    form_ = FormDiagram.from_lines(lines + symmetry, delete_boundary_face=False)
    form_.update_default_vertex_attributes({'is_roller': False})
    form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    form_.attributes['loadpath'] = 0

    gkey_key = form_.gkey_key()

    for i in pins:
        form_.set_vertex_attribute(gkey_key[i], 'is_fixed', value=True)

    for i in symmetry:
        a, b = i
        u = gkey_key[geometric_key(a)]
        v = gkey_key[geometric_key(b)]
        form_.set_edge_attribute((u, v), name='is_symmetry', value=True)

    # Loads

    pz = 0
    for key in form_.vertices():
        gkey = geometric_key(form_.vertex_coordinates(key))
        form_.vertex[key]['pz'] = loads[gkey]
        pz += form_.vertex[key]['pz']
    print('Form: ' + ST)
    print('Total load: {0}'.format(pz))

    plot_form(form_).show()
    form_.to_json(radical + '_sym.json')
