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

from copy import deepcopy

import math

from compas_thrust.plotters.plotters import plot_form

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'not_sym_load',
    'fill_load',
    'set_cross_vault_loads',
    'set_pavillion_vault_loads',
    'set_oct_vault_loads',
    'set_dome_loads',
]

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

def fill_load(form, t = 0.5, ro = 1.0, scale = 100):

    z = [form.vertex_coordinates(key)[2] for key in form.vertices()]
    zmax = max(z)
    print('Max Z of form is: {0:.1f}'.format(zmax))
    pzt = 0

    form0 = deepcopy(form)
    for key in form0.vertices():
        form0.set_vertex_attribute(key, 'z', value = 0.0)
    z0 = [form0.vertex_coordinates(key)[2] for key in form0.vertices()]

    for key in form.vertices():
        _, _, zi = form.vertex_coordinates(key)
        ai = form0.vertex_area(key=key)
        hi = zmax - zi
        pzi = (hi + t) + ro * ai
        form.set_vertex_attribute(key, 'pz', value = pzi)
        pzt += pzi

    print('Total Load (filling) applied equals: {0:.2f}'.format(pzt))

    if scale:
        scl = scale/pzt
        for key in form.vertices():
            form.vertex[key]['pz'] *= scl

        print('Load scaled to a total of: {0:.2f}'.format(pzt*scl))


    return form

def set_cross_vault_loads(form, xy_span = [[0.0,10.0],[0.0,10.0]], thickness = None, tol = 0.00, set_heights = False):

    # WIP

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi <= y1/x1 * xi + tol and yi <= y1 - xi + tol: #Q1
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi >= y1/x1 * xi - tol and yi >= y1 - xi - tol: #Q3
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi <= y1/x1 * xi + tol and yi >= y1 - xi - tol: #Q2
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        elif yi >= y1/x1 * xi - tol and yi <= y1 - xi + tol: #Q4
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    return form

def set_pavillion_vault_loads(form, xy_span = [[0.0,10.0],[0.0,10.0]], thickness = None, tol = 0.00, set_heights = False):

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi <= y1/x1 * xi + tol and yi <= y1 - xi + tol: #Q1
            try:
                z = math.sqrt((ry)**2 - (yi-ry)**2)
            except:
                z = 0
        elif yi >= y1/x1 * xi - tol and yi >= y1 - xi - tol: #Q3
            try:
                z = math.sqrt((ry)**2 - (yi-ry)**2)
            except:
                z = 0
        elif yi <= y1/x1 * xi + tol and yi >= y1 - xi - tol: #Q2
            try:
                z = math.sqrt((rx)**2 - (xi-rx)**2)
            except:
                z = 0
        elif yi >= y1/x1 * xi - tol and yi <= y1 - xi + tol: #Q4
            try:
                z = math.sqrt((rx)**2 - (xi-rx)**2)
            except:
                z = 0
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)
        if set_heights:
            form.set_vertex_attribute(key,'z',value=round(z,2))

    return form

def set_oct_vault_loads(form, xy_span = [[0.0,10.0],[0.0,10.0]], thickness = None, tol = 0.01):

    # WIP
    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0])/2.0

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if yi < y1/x1*xi + tol and yi < y1 - x1 + tol: #Q1
            z = math.sqrt= ((rx)**2 - (xi-rx)**2)
        elif yi > y1/x1*xi - tol and yi > y1 - x1 - tol: #Q3
            z = math.sqrt((rx)**2 - (xi-rx)**2)
        elif yi < y1/x1*xi + tol and yi > y1 - x1 - tol: #Q2
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        elif yi > y1/x1*xi - tol and yi < y1 - x1 + tol: #Q4
            z = math.sqrt((ry)**2 - (yi-ry)**2)
        else:
            print('Vertex {0} did not belong to any Q. (x,y) = ({1},{2})'.format(key,xi,yi))
            z = 0.0
        form.set_vertex_attribute(key,'target',value=z)

    return form

def set_dome_loads(form, center = [0.0,0.0], radius = 10.0, thickness = None, tol = 0.00, set_heights = False, scale =100):

    x0 = center[0]
    y0 = center[1]

    loads = deepcopy(form)

    for key in loads.vertices():
        xi, yi, _ = loads.vertex_coordinates(key)
        z2 = + radius**2 - (xi - x0)**2 - (yi - y0)**2
        if -0.01 <= z2 <= 0.0:
            z2 = 0.0
        try: 
            z = math.sqrt(z2)
        except:
            print(xi,yi)
            z=0
        loads.set_vertex_attribute(key,'z',value=z)
    
    pzt = 0.0
    for key in form.vertices():
        pz = loads.vertex_area(key)
        form.set_vertex_attribute(key,'pz',value=pz)
        pzt += pz

    print('Total_Load: {0:.2f}'.format(pzt))


    if scale:
        scl = scale/pzt
        pzt = 0.0
        for key in form.vertices():
            form.vertex[key]['pz'] *= scl
            pzt += form.vertex[key]['pz']

    print('After Scaling Load: {0:.2f}'.format(pzt))

    return form
