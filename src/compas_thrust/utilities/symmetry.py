from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas.utilities import geometric_key    
from compas.utilities import reverse_geometric_key

from compas.geometry.transformations.transformations import mirror_point_line
from compas.geometry.transformations.transformations import rotate_points
from compas.geometry import closest_point_on_line
from compas.geometry import midpoint_line_xy
from compas.geometry import matrix_from_axis_and_angle

from compas.geometry.distance import distance_point_point_xy
from numpy import argmin
from numpy import sqrt
from compas_thrust.algorithms.equilibrium import z_from_form
from compas_thrust.diagrams.form import overview_forces

from compas.datastructures import Mesh
from compas_plotters import MeshPlotter

import math

from compas_thrust.plotters.plotters import plot_form

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'replicate',
    'replicate2',
    'create_sym',
    'create_sym2',
    'fix_boundaries_sym',
    'fix_boundaries_complete',
    'fix_mid_sym',
    'fix_mid_complete',
    'not_sym_load',
]

def replicate(form,file, plot=False):

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
            gkey_appx = geometric_key(real_points[point_i])
            form_.set_edge_attribute((u, v), name='q', value = q_i[gkey_appx])

        if gkey_mid in mid_ind:
            form_.set_edge_attribute((u, v), name='is_ind', value = True)
        else:
            form_.set_edge_attribute((u, v), name='is_ind', value = False)


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
    if plot:
        overview_forces(form_)

    if plot:
        plot_form(form_, radius=0.05).show()

    return form_

def replicate2(form, file, plot=None):

    """ Copy 1/2 form-diagram in a complete structure.

    Parameters
    ----------
    form : obj
        FormDiagram 1/8.

    form : json
        Location of the form-diagram.

    plot : bool
        Plot the Initial and Final configuration.

    Returns
    -------
    form : obj
        complete FormDiagram.

    """

    form_ = FormDiagram.from_json(file)

    if plot:
        plot_form(form_, radius=0.05).show()

    q_i = {}
    mid_ind = []
    tol = 0.001

    for u, v in form.edges():
        if form.get_edge_attribute((u,v),'is_symmetry') == False:
            sp, ep = form.vertex_coordinates(u), form.vertex_coordinates(v)
            line_edge = [sp, ep]
            mp = midpoint_line_xy(line_edge)
            gkey_mid = geometric_key(mp)
            q_i[gkey_mid] = form.get_edge_attribute((u,v),'q')

            line_sym = ([0.0,0.0,0.0],[0.0,10.0,0.0])

            mirror_diag = mirror_point_line(mp,line_sym)
            gkey_mirror = geometric_key(mirror_diag)

            q_i[gkey_mirror] = q_i[gkey_mid]
            
            if -tol <= sp[0] <= tol and -tol <= ep[0] <= tol:
                q_i[gkey_mirror] *= 2

        if form.get_edge_attribute((u,v),'is_ind') == True:
            sp, ep = form.vertex_coordinates(u), form.vertex_coordinates(v)
            line_edge = [sp, ep]
            mp = midpoint_line_xy(line_edge)
            gkey_mid = geometric_key(mp)
            mid_ind.append(gkey_mid)

    real_points = []
    for key in q_i:
        real_points.append(reverse_geometric_key(key))

    for u, v in form_.edges():
        sp, ep = form_.vertex_coordinates(u), form_.vertex_coordinates(v)
        line_edge = [sp, ep]
        mp = midpoint_line_xy(line_edge)
        gkey_mid = geometric_key(mp)
        try:
            form_.set_edge_attribute((u, v), name='q', value = q_i[gkey_mid])
        except:
            dist = []
            for point in real_points:
                dist.append(distance_point_point_xy(point,mp))
            point_i = argmin(dist)
            gkey_appx = geometric_key(real_points[point_i])
            form_.set_edge_attribute((u, v), name='q', value = q_i[gkey_appx])

        if gkey_mid in mid_ind:
            form_.set_edge_attribute((u, v), name='is_ind', value = True)
        else:
            form_.set_edge_attribute((u, v), name='is_ind', value = False)

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
    overview_forces(form_)

    if plot:
        plot_form(form_, radius=0.05).show()

    return form_

def create_sym(form, keep_q = True):

    plot_form(form).show()

    lines = []
    symmetry = []
    pins = []
    loads = {}
    qs = {}
    target = {}
    lb = {}
    ub = {}
    tol = 0.001

    pz = 0
    for key in form.vertices():
        pz += form.vertex[key]['pz']
    print('Load Before Sym: {0}'.format(pz))

    for u,v in form.edges():
        u_coord = form.vertex_coordinates(u)[:2] + [0]
        v_coord = form.vertex_coordinates(v)[:2] + [0]
        mid_gkey = geometric_key(form.edge_midpoint(u,v)[:2] + [0])
        qs[mid_gkey] = form.get_edge_attribute((u,v), 'q')
        if u_coord[1] >= 10 - u_coord[0] - tol and u_coord[0] <= 5.0 + tol and v_coord[1] >= 10 - v_coord[0] - tol and v_coord[0] <= 5.0 + tol :
            lines.append([u_coord,v_coord])
        if 5.0 - tol <= u_coord[0] <= 5.0 + tol and 5.0 - tol <= v_coord[0] <= 5.0 + tol: # Line Vertical at x = 5.0 and y > 5.0
            qs[mid_gkey] *= 0.5
        if 10 - u_coord[0] - tol <= u_coord[1] <= 10 - u_coord[0] + tol and 10 - v_coord[0] - tol <= v_coord[1] <= 10 - v_coord[0] + tol: # Line Diagonal at y = x 
            qs[mid_gkey] *= 0.5

    for key in form.vertices():
        coord = form.vertex_coordinates(key)[:2] + [0]
        gkey = geometric_key(coord)
        loads[gkey] = form.get_vertex_attribute(key, 'pz')
        target[gkey] = form.get_vertex_attribute(key, 'target')
        lb[gkey] = form.get_vertex_attribute(key, 'lb')
        ub[gkey] = form.get_vertex_attribute(key, 'ub')

        if 5.0 - tol <= coord[0] <= 5.0 + tol and coord[1] > 5.0 - tol: # Line Vertical at x = 5.0 and y > 5.0
            loads[gkey] = 0.5 * loads[gkey]
            if form.get_vertex_attribute(key, 'is_fixed') == False:
                ext_pt = [coord[0]+1.0, coord[1], 0.0]
                symmetry.append([coord,ext_pt])
                pins.append(geometric_key(ext_pt))
                loads[geometric_key(ext_pt)] = 0.0

        if 10 - coord[0] - tol <= coord[1] <= 10 - coord[0] + tol and 10.0 - tol > coord[1] > 5.0 + tol: # Line Diagonal at y = x - 10 without central and corner
            loads[gkey] = 0.5 * loads[gkey]
            ext_pt = [coord[0]-sqrt(2)/2,coord[1]-sqrt(2)/2,0.0]
            symmetry.append([coord,ext_pt])
            pins.append(geometric_key(ext_pt))
            loads[geometric_key(ext_pt)] = 0.0

        if 5.0 - tol <= coord[0] <= 5.0 + tol and 5.0 - tol <= coord[1] <= 5.0 + tol: # Central Point
            loads[gkey] = 0.25 * loads[gkey]
            ext_pt = [coord[0],coord[1]-1.0,0.0]
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

    plot_form(form_).show()

    for i in pins:
        form_.set_vertex_attribute(gkey_key[i], 'is_fixed', value=True)

    for i in symmetry:
        a, b = i
        u = gkey_key[geometric_key(a)]
        v = gkey_key[geometric_key(b)]
        form_.set_edge_attribute((u, v), name='is_symmetry', value=True)

    # Loads and Qs

    pz = 0
    for key in form_.vertices():
        gkey = geometric_key(form_.vertex_coordinates(key))
        form_.vertex[key]['pz'] = loads[gkey]
        pz += form_.vertex[key]['pz']
        try:
            form_.set_vertex_attribute(key, name = 'target', value = target[gkey])
            form_.set_vertex_attribute(key, name = 'lb', value = lb[gkey])
            form_.set_vertex_attribute(key, name = 'ub', value = ub[gkey])
        except:
            pass
    print('Load After Sym: {0}'.format(pz))

    for u, v in form_.edges_where({'is_symmetry': False}):
        qi = qs[geometric_key(form_.edge_midpoint(u,v))]
        form_.set_edge_attribute((u, v), name = 'q', value = qi)

    plot_form(form_).show()

    return form_

def create_sym2(form, keep_q = True):

    lines = []
    symmetry = []
    pins = []
    loads = {}
    qs = {}
    tol = 0.001

    for key, attr in form.vertices(True):
        attr['z'] = 0.0

    for u,v in form.edges():
        u_coord = form.vertex_coordinates(u)
        v_coord = form.vertex_coordinates(v)
        qs[geometric_key(form.edge_midpoint(u,v))] = form.get_edge_attribute((u,v), 'q')
        if u_coord[0] <= 0.0 + tol and v_coord[0] <= 0.0 + tol: # Line Vertical at x = 0.0
            qs[geometric_key(form.edge_midpoint(u,v))] *= 0.5 # isn't it changing all q's for half? VERIFY
            lines.append([u_coord,v_coord])

    for key in form.vertices():
        coord = form.vertex_coordinates(key)
        gkey = geometric_key(coord)
        loads[gkey] = form.get_vertex_attribute(key, 'pz')

        if 0.0 - tol <= coord[0] <= 0.0 + tol: # Line Vertical at x = 0.0
            loads[gkey] = 0.5 * loads[gkey]
            if form.get_vertex_attribute(key, 'is_fixed') == False:
                ext_pt = [coord[0]+1.0, coord[1], coord[2]]
                symmetry.append([coord,ext_pt])
                pins.append(geometric_key(ext_pt))
                loads[geometric_key(ext_pt)] = 0.0

        if coord[0] <= 0.0 + tol and form.get_vertex_attribute(key, 'is_fixed') == True: # Pins To the Left
            pins.append(geometric_key(coord))


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

    for u, v in form_.edges_where({'is_symmetry': False}):
        qi = qs[geometric_key(form_.edge_midpoint(u,v))]
        form_.set_edge_attribute((u, v), name = 'q', value = qi)

    plot_form(form_).show()

    # Loads

    pz = 0
    for key in form_.vertices():
        gkey = geometric_key(form_.vertex_coordinates(key))
        form_.vertex[key]['pz'] = loads[gkey]
        pz += form_.vertex[key]['pz']
    print('Total load: {0}'.format(pz))

    return form_

def fix_boundaries_sym(form, plot = False):

    tol = 0.001

    for key in form.vertices():
        _, y, _ = form.vertex_coordinates(key)
        if y > 10.0 - tol and y < 10.0 + tol:
            form.set_vertex_attribute(key, 'is_fixed', True)
            form.set_vertex_attribute(key, 'z', 0.0)

    if plot:
        plot_form(form).show()

    return form

def fix_boundaries_complete(form, plot = False):

    for key in form.vertices_on_boundary():
        form.set_vertex_attribute(key, 'is_fixed', True)
        form.set_vertex_attribute(key, 'z', 0.0)

    if plot:
        plot_form(form).show()

    return form

def fix_mid_sym(form, plot = False):

    tol = 0.001

    for key in form.vertices():
        x, y, _ = form.vertex_coordinates(key)
        if y > 10.0 - tol and y < 10.0 + tol and x > 5.0 - tol and x < 5.0 + tol:
            form.set_vertex_attribute(key, 'is_fixed', True)
            form.set_vertex_attribute(key, 'z', 0.0)

    if plot:
        plot_form(form, show_q=False, fix_width=True).show()

    return form

def fix_mid_complete(form, plot = False):

    tol = 0.001

    for key in form.vertices_on_boundary():
        x, y, _ = form.vertex_coordinates(key)
        if ( (y > 10.0 - tol and y < 10.0 + tol) or  (y > 0.0 - tol and y < 0.0 + tol) ) and (x > 5.0 - tol and x < 5.0 + tol):
            form.set_vertex_attribute(key, 'is_fixed', True)
            form.set_vertex_attribute(key, 'z', 0.0)
        if ( (x > 10.0 - tol and x < 10.0 + tol) or  (x > 0.0 - tol and x < 0.0 + tol) ) and (y > 5.0 - tol and y < 5.0 + tol):
            form.set_vertex_attribute(key, 'is_fixed', True)
            form.set_vertex_attribute(key, 'z', 0.0)

    if plot:
        plot_form(form, show_q=False, fix_width=True).show()

    return form

def not_sym_load(form, x0 = 0, x1 = 5.0, magnitude = 2.0):

    tol = 0.01

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        if x > x0 - tol and x < x1 + tol:
            pz0 = form.get_vertex_attribute(key, 'pz')
            if x > x1 - tol:
                form.set_vertex_attribute(key, 'pz', value = ((magnitude-1)/2 +1) * pz0)
            else:
                form.set_vertex_attribute(key, 'pz', value = magnitude * pz0)
    
    return form


