"""
Script to compare the time consumption for finding the independents
"""

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.problems import initialise_form
from compas_tno.problems import initialise_problem
from compas_tno.algorithms import find_independents
from compas.numerical import rref_sympy
from compas.numerical import nonpivots
from compas.numerical import connectivity_matrix
from compas.numerical import equilibrium_matrix
from numpy.linalg import matrix_rank
import time

# type_formdiagrams = ['cross_fd', 'fan_fd']
# discretisations = [[10, 12, 14, 16, 18, 20], [10, 12, 14, 16, 18, 20]]
# type_formdiagrams = ['radial_fd']
discretisations = [[[24, 12], [24, 16], [24, 20], [24, 24]]]
# discretisations = [[[4, 12], [4, 16]]]

type_formdiagrams = ['cross_fd']
discretisations = [[10, 12, 14]]

sols = {}

for i in range(len(type_formdiagrams)):
    type_formdiagram = type_formdiagrams[i]
    sols[type_formdiagram] = {}

    for j in range(len(discretisations[i])):
        discretisation = discretisations[i][j]
        sols[type_formdiagram][str(discretisation)] = {}
        print(discretisation)

        span = 10.0

        data_diagram = {
            'type': type_formdiagram,
            'xy_span': [[0, span], [0, span]],
            'discretisation': discretisation,
            'fix': 'corners',
        }

        if type_formdiagram == 'radial_fd':
            radius = 5.0
            data_diagram['discretisation'] = [discretisation[0], discretisation[1]]
            data_diagram['center'] = [5.0, 5.0]
            data_diagram['radius'] = radius

        form = FormDiagram.from_library(data_diagram)

        type_structure = 'crossvault'
        thk = 0.50
        data_shape = {
            'type': type_structure,
            'thk': thk,
            'discretisation': discretisation,
            'xy_span': [[0, span], [0, span]],
            't': 0.0,
        }

        if type_formdiagram == 'radial_fd':
            data_shape['type'] = 'dome'
            radius = 5.0
            data_shape['discretisation'] = [discretisation[0], discretisation[1]]
            data_shape['center'] = [5.0, 5.0]
            data_shape['radius'] = radius

        shape = Shape.from_library(data_shape)

        # from compas_tno.viewers import view_shapes
        # import matplotlib
        # import matplotlib.pyplot as plt
        # matplotlib.use('TkAgg')
        # import os
        # os.environ['QT_MAC_WANTS_LAYER'] = '1'
        # view_shapes(shape).show()

        form.selfweight_from_shape(shape)

        # initialise_form(form, printout=True)
        # initialise_form(form, printout=True)

        vertex_index = form.vertex_index()

        xy = form.vertices_attributes('xy')
        fixed = [vertex_index[vertex] for vertex in form.fixed()]
        free = list(set(range(form.number_of_vertices())) - set(fixed))
        edges = [(vertex_index[u], vertex_index[v]) for u, v in form.edges()]
        C = connectivity_matrix(edges)
        E = equilibrium_matrix(C, xy, free)

        print('-'*10, 'Problem type:', type_formdiagram, 'Discretisation:', discretisation)

        print('Shape of E:', E.shape)
        mr = matrix_rank(E)
        print('Rank E:', mr)
        print('Expected ind. count:', mr - E.shape[1])

        start_time = time.time()
        ind = find_independents(E)
        elapsed_time = time.time() - start_time
        print('Elapsed Time with TNO: {0:.1f} sec'.format(elapsed_time))
        ind.sort()
        print('TNO found inds:', len(ind))
        print('TNO ind:', ind)
        sols[type_formdiagram][str(discretisation)]['tno'] = round(elapsed_time, 2)

        start_time = time.time()
        ind = nonpivots(rref_sympy(E))
        elapsed_time = time.time() - start_time
        print('-')
        print('Elapsed Time with COMPAS: {0:.1f} sec'.format(elapsed_time))
        sols[type_formdiagram][str(discretisation)]['ags'] = round(elapsed_time, 2)
        ind.sort()
        print('AGS found inds:', len(ind))
        print('AGS ind:', ind)

    print(sols)

print('******** Resume of Solutions')
for key in sols:
    print(key)
    for discr in sols[key]:
        print(discr)
        for type_algo in sols[key][discr]:
            print(type_algo, ',', sols[key][discr][type_algo])
