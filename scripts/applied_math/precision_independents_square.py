import os
os.environ["USE_PROPACK"]="1"

from scipy.sparse.linalg import svds
from numpy.linalg import svd
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt
import numpy as np
from numpy import asarray

from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas.datastructures import Mesh
from compas_tno.problems import initialise_form
from compas_tno.utilities import slide_diagram
from compas_tno.algorithms import check_independents
import numpy as np

discretisation = 10
form = FormDiagram.create_cross_form(discretisation=discretisation)

for lambd in [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:

    form = FormDiagram.create_parametric_form(discretisation=discretisation, lambd=lambd)

    # # ---------- EXAMPLE 1

    M = initialise_form(form, find_inds=True, printout=True)
    k = len(M.ind)
    print(M.ind)

    n, m = M.E.shape
    print(M.E.shape)

    plotter = TNOPlotter(form)
    plotter.draw_form_independents()
    plotter.show()

    _, s, _ = svd(asarray(M.E))
    # _, s, _ = svds(M.E, k=min(M.E.shape), solver='propack')
    print('max/min singular vectors E', max(s), min(s), len(s))

    mn = max(m - n, n - m)
    zs = k - mn
    print('ZS:', zs)
    print(k - mn)
    print('Non zero SVs:', s[:len(s)-zs][:6], '...')
    print('Zero SVs:', s[-zs:])

    check = check_independents(M)
    print('Check independents:', check)

    if zs > 0:
        last_zero = s[len(s)-zs]
        first_non_zero = s[len(s)-zs - 1]

        print('First Non Zero:', first_non_zero)
        print('Last Zero SV:', last_zero)

        len_zero = len(s[s<1.0])
        lin_x = len_zero - zs

        porcentage_key = (first_non_zero - last_zero)/first_non_zero
        print('percentage key is:', porcentage_key)

        fig, ax = plt.subplots()
        ax.plot(s[s<1.0])
        ax.plot([lin_x, lin_x], [0, 1.0], color='black')
        plt.show()

    else:
        fig, ax = plt.subplots()
        ax.plot(s[s<1.0])
        plt.show()

# # ---------- EXAMPLE 3

import compas_tno

for type_mesh in ['G2']:  # ['A2', 'B2', 'C2', 'D2', 'E2']:

    print('---- Mesh ', type_mesh)

    path = '/Users/mricardo/compas_dev/me/pattern/singular/dome/mesh-' + type_mesh + '.json'
    path = '/Users/mricardo/compas_dev/me/max_load/dome/apex/dome/mesh-' + type_mesh + '.json'

    # ad = compas_tno.get('form.json')
    # form = FormDiagram.from_json(ad)
    mesh = Mesh.from_json(path)
    form = FormDiagram.from_mesh(mesh)
    form.delete_boundary_edges()
    form.set_boundary_supports()

    plotter = TNOPlotter(form)
    plotter.draw_form()
    plotter.show()

    M = initialise_form(form, find_inds=True, printout=True)
    k = len(M.ind)
    print('Independents Robin', k)

    n, m = M.E.shape
    print(M.E.shape)
    print('rank E:', matrix_rank(M.E))
    _, s, _ = svds(M.E, k=min(M.E.shape), solver='propack')
    print('max/min singular vectors E', max(s), min(s), len(s))

    zs = k - (m - n)
    print(k - (m - n))
    print('Non zero SVs:', s[zs:][:6], '...')
    print('Zero SVs:', s[:zs])

    if zs > 0:

        last_zero = s[:zs][-1]
        first_non_zero = s[zs:][0]
        print('First Non Zero:', first_non_zero)
        print('Last Zero SV:', last_zero)

        porcentage_key = (first_non_zero - last_zero)/first_non_zero
        print('percentage key is:', porcentage_key)

        fig, ax = plt.subplots()
        # ax.plot(s[::-1])
        ax.plot(s[s<1.0][::-1])
        plt.show()

    else:
        fig, ax = plt.subplots()
        ax.plot(s[s<1.0][::-1])
        plt.show()

    plotter = TNOPlotter(form)
    plotter.draw_form_independents()
    plotter.show()

    check = check_independents(M)

    print('Check independents:', check)

    # path_inds = '/Users/mricardo/compas_dev/me/pattern/singular/dome/mesh-' + type_mesh + '-with_inds.json'
    # form.to_json(path_inds)
