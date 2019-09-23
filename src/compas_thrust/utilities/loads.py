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

