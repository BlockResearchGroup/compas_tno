import os
os.environ["USE_PROPACK"]="1"

from turtle import shape
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_plotters import Plotter
from compas.datastructures import Mesh
from compas_tno.problems import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.algorithms import xyz_from_xopt
from numpy import zeros
from numpy import array
import compas_tno
import os
from copy import deepcopy


def check_ind_sensitivity(form, M):

    ind = M.ind
    ind_good = []
    inds = {}
    text = {}
    i = 0
    for edge in form.edges_where({'_is_edge': True}):
        if form.edge_attribute(edge, 'is_ind'):
            text[edge] = str(i)
            inds[i] = edge
            i = i+1

    z0 = form.vertices_attributes('z')
    z0 = array(z0).reshape(-1, 1)

    q_init = [M.q[i] for i in range(len(M.q))]

    delta = -10

    inds_problem = []
    inds_ok = []
    inds_problem_print = {}
    inds_ok_print = {}

    for i_ind in range(len(M.ind)):
        print('\n', '-'*10)
        print('i=', i_ind)
        M.q = array(q_init).reshape(-1, 1)
        print('max min q - before:', max(M.q.flatten()), min(M.q.flatten()))

        edge = inds[i_ind]
        q0 = M.q[ind[i_ind]]
        print('q0 = ', q0)

        M.q[ind[i_ind]] += delta
        print('new q0 = ', M.q[ind[i_ind]])

        print('variables: ', M.q[ind].flatten())
        M = xyz_from_xopt(M.q[ind], M)

        for i, vertex in enumerate(form.vertices()):
            # x, y, z = M.X[i]
            form.vertex_attributes(vertex, 'xyz', M.X[i])

        for i, edge in enumerate(form.edges_where({'_is_edge': True})):
            form.edge_attribute(edge, 'q', M.q[i][0])

        print('max min q:', max(M.q.flatten()), min(M.q.flatten()))

        dz = (M.X[:, 2].reshape(-1, 1) - z0) ** 2
        dq = (M.q - array(q_init).reshape(-1, 1)) ** 2

        print('sum dz:', sum(dz))
        print('sum dq:', sum(dq))

        if sum(dq)[0] > 1e+10:
            inds_problem.append(i_ind)
            inds_problem_print[inds[i_ind]] = i_ind
        else:
            inds_ok_print[inds[i_ind]] = i_ind
            inds_ok.append(i_ind)
            ind_good.append(ind[i_ind])

    print('inds ok', inds_ok, len(inds_ok))
    print('inds not ok', inds_problem, len(inds_problem))

    print('Showing inds OK')

    plotter = TNOPlotter(form, shape=shape)
    plotter.draw_supports()
    plotter.draw_form_independents()
    plotter.formartist.draw_edgelabels(text=inds_ok_print)
    # plotter.draw_form()
    plotter.show()

    print('Showing inds not OK')

    plotter = TNOPlotter(form, shape=shape)
    plotter.draw_supports()
    plotter.draw_form_independents()
    plotter.formartist.draw_edgelabels(text=inds_problem_print)
    # plotter.draw_form()
    plotter.show()

    return inds_ok


# ad_analysis = compas_tno.get('analysis.json')
# ad_analysis = compas_tno.get('analysis-domesimple.json')
# ad_analysis = compas_tno.get('analysis-domeC.json')
# ad_analysis = compas_tno.get('analysis-domeD.json')
# ad_analysis = compas_tno.get('analysis-domeA2.json')
# ad_analysis = compas_tno.get('analysis-domeB2tna.json')
# ad_analysis = compas_tno.get('analysis-cross.json')
# ad_analysis = compas_tno.get('analysis-cross14.json')
# ad_analysis = compas_tno.get('analysis-cross-sag.json')

# robin
ad_analysis = compas_tno.get('analysis-domeA2.json')
# ad_analysis = compas_tno.get('analysis-domeB2.json')
# ad_analysis = compas_tno.get('analysis-domeC2.json')
# ad_analysis = compas_tno.get('analysis-domeD2.json')
# ad_analysis = compas_tno.get('analysis-domeE2.json')

analysis = Analysis.from_json(ad_analysis)

analysis.optimiser.settings['variables'] = ['q']
analysis.optimiser.settings['starting_point'] = 'current'
analysis.optimiser.settings['plot'] = False

analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

form = analysis.form
shape = analysis.shape
print(shape)

M = analysis.optimiser.M

# viewer = Viewer(form, shape)
# viewer.draw_thrust()
# viewer.draw_shape()
# viewer.draw_force()
# viewer.show()

ind_good = check_ind_sensitivity(form, M)
ind_bad = list(set(range(len(M.ind))) - set(ind_good))

from numpy.linalg import cond
condnumber = cond(M.B)
print('Cond number B is:', condnumber)
condnumber2 = cond(M.B[:, ind_good])
print('Cond number B[ind_good] 2 is:', condnumber2)
nl, nc = M.B.shape
print(M.B.shape)
print(M.B[:, ind_good].shape)

# column_bad = M.B[:, ind_bad[0]].flatten()
# column_good = M.B[:, ind_good[0]].flatten()

# for i in range(nc):
#     print('\n-----------------\ncolumn i:', i)
#     column = M.B[:, i].flatten()
#     if i in ind_good:
#         text = '(good)'
#     else:
#         text = '(bad)'
#     print(text + 'max, min:', max(column), min(column))

# print('column bad | max | min | vector:', max(column_bad), min(column_bad), column_bad)
# print('column good | max | min | vector:', max(column_good), min(column_good), column_good)


from scipy.sparse.linalg import svds
from scipy.linalg import diagsvd
from numpy import dot
from numpy import diag  # maybe diags?
from numpy.linalg import matrix_rank  # maybe diags?

m, n = M.E.shape
U, s, Vh = svds(M.E, k=min(m, n), solver='propack')

print('Cond number E is:', condnumber)
rank = matrix_rank(M.E)
print('Rank E is:', rank)
print('Shape E is:', M.E.shape)
print('Shape S is:', s.shape)
print('Shape U is:', U.shape)
print('Shape Vh is:', Vh.shape)
print('Number of inds is:', len(M.ind))


# Edinv = array(M.Edinv)
# print('Edinv.shape =', Edinv.shape)
# print(type(Edinv), type(M.B))
# condnumber = cond(Edinv)
# print('Edinv Cond number is:', condnumber)
# condnumber2 = cond(Edinv[:, ind_good])
# print('Edinv Cond number 2 is:', condnumber2)

disc = len(M.ind) - (n - m)
disc_val = s[:disc]
print('E sing. values:', s)
print('max/min singular vectors E', max(s), min(s), len(s))
print('Discated singular vectors E', disc, disc_val)
print('Maximum Discated Singular Value:', max(disc_val), '< 1e-4?')
print('First considered Singular Value:', s[disc], 'much bigger?')

# m, n = M.E[:, M.dep].shape
# U, s, Vh = svds(M.E[:, M.dep], k=min(m, n), solver='propack')
# print('\nEd sing. values:', s)
# print('max/min singular vectors Ed', max(s), min(s), len(s))

# m, n = M.Edinv.shape
# U, s, Vh = svds(M.Edinv, k=min(m, n), solver='propack')
# print('\nEdinv sing. values:', s)
# print('max/min singular vectors Edinv', max(s), min(s), len(s))


# # Solving the system (useless)
# c = dot(U.T, M.ph)
# w = dot(diag(1/s), c)
# x = dot(Vh.T, w)
# print(max(x), min(x))

# # Geometry parameters

# radius = 5.0
# thk = 0.5
# discretisation = [8, 10]
# discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

# # Parameters Optimisations

# obj = 'max_load'
# solver = 'IPOPT'
# constraints = ['funicular', 'envelope', 'reac_bounds']
# variables = ['q', 'zb', 'lambdv']  # lambdv
# features = ['fixed']
# starting_point = 'loadpath'
# make_video = True
# autodiff = False

# # Create shape/diagram

# dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.0)

# ad_form = compas_tno.get('form-lp.json')
# form = FormDiagram.from_json(ad_form)

# optimiser = Optimiser()
# optimiser.settings['objective'] = obj
# optimiser.settings['solver'] = solver
# optimiser.settings['constraints'] = constraints
# optimiser.settings['variables'] = variables
# optimiser.settings['features'] = features
# optimiser.settings['starting_point'] = starting_point
# optimiser.settings['derivative_test'] = False
# optimiser.settings['printout'] = True
# optimiser.settings['plot'] = True
# optimiser.settings['save_iterations'] = make_video
# optimiser.settings['autodiff'] = autodiff

# max_load_mult = 600.0
# n = form.number_of_vertices()
# pzv = zeros((n, 1))
# # pzv[0] = -1.0
# pzv[48] = -1.0
# # pzv[37] = -1.0

# optimiser.settings['max_lambd'] = max_load_mult
# optimiser.settings['load_direction'] = pzv

# # Create analysis

# print(type(optimiser.settings['load_direction']))

# analysis = Analysis.from_elements(dome, form, optimiser)

# ad_shape = compas_tno.get('shape.json')
# dome.to_json(ad_shape)
# print('File saved @', ad_shape)

# ad_optimiser = compas_tno.get('optimiser.json')
# optimiser.to_json(ad_optimiser)
# print('File saved @', ad_optimiser)

# ad_analysis = compas_tno.get('analysis.json')
# analysis.to_json(ad_analysis)
# print('File saved @', ad_analysis)

# form = FormDiagram.from_json(ad_form)
# shape = Shape.from_json(ad_shape)
# optimiser = Optimiser.from_json(ad_optimiser)
# analysis = Analysis.from_json(ad_analysis)

# print(type(optimiser.settings['load_direction']))

# print(optimiser.M)

# # analysis.apply_selfweight()
# # analysis.apply_envelope()
# # analysis.apply_reaction_bounds()
# # analysis.set_up_optimiser()
# # analysis.run()


