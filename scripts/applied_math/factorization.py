from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from numpy import identity, zeros, ones
from numpy.linalg import matrix_rank
from numpy.linalg import pinv

import scipy

# **************************************************
#
# SCRIPT TO CHECK IF IT IS POSSIBLE TO SPEED UP INDEP. FINDING WITH LU/QR FACTORIZATIONS
#
# **************************************************

# Geometry parameters

radius = 5.0
thk = 0.5
discretisation = [8, 10]
discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

# Parameters Optimisations

obj = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'loadpath'

# Create shape/diagram

dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.5)

form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)

# Create optimiser

optimiser = Optimiser()
optimiser.settings['objective'] = obj
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['starting_point'] = starting_point
optimiser.settings['derivative_test'] = True
optimiser.settings['printout'] = True

# Create analysis

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

# import json
# from numpy import array
# from numpy import zeros

# import cvxpy as cp
# from cvxpy import diag
# from cvxpy import matrix_frac
# from cvxpy import Minimize
# from cvxpy import Problem

# path = '/Users/mricardo/compas_dev/compas_tno/data/data.json'

# with open(path, 'r') as json_file:
#     data = json.load(json_file)

# q = array(data['q'])
# E = array(data['E'])
# C = array(data['C'])
# Ci = array(data['Ci'])
# Cit = array(data['Cit'])
# Cb = array(data['Cb'])
# pz = array(data['pz'])
# qmin = array(data['qmin'])
# qmax = array(data['qmax'])
# free = data['free']
# fixed = data['fixed']
# x = array(data['x'])
# y = array(data['y'])
# m = data['m']

M = optimiser.M
m = M.m

E = M.E
ph = M.ph

# i = 0
# ind_ = []
# tol = 1e-10

# Q, R = scipy.linalg.qr(E)

# for ncol in range(R.shape[1]):
#     el = R[i][ncol]
#     print(i, ncol, el)
#     if abs(el) < tol:
#         ind_.append(ncol)
#         print('Added i ncol', ncol)
#     else:
#         i = i + 1

# print(ind_)
# print(M.ind)

i = 0
ind_ = []
tol = 1e-10

P, L, U = scipy.linalg.lu(E)

print('U shape', U.shape, matrix_rank(U))

for ncol in range(U.shape[1]):
    el = U[i][ncol]
    print(i, ncol, el)
    if abs(el) < tol:
        print('Found not-pivot at (i, ncol)', i, ncol)
        for nlin in ind_:
            if nlin < U.shape[0]:
                el_ = U[nlin][ncol]
                print('Checking for additional pivot at (nlin, ncol)', nlin, ncol, el_)
        ind_.append(ncol)

    if i < U.shape[0] - 1:
        i += 1


print(ind_)
dep_ = list(set(range(U.shape[1])) - set(ind_))
print(M.ind)

Uind = U[:, M.ind]
Udep = U[:, M.dep]
Uind_ = U[:, ind_]
Udep_ = U[:, dep_]

print(Uind.shape, matrix_rank(Uind))
print(Uind_.shape, matrix_rank(Uind_))

print(Udep.shape, matrix_rank(Udep))
print(Udep_.shape, matrix_rank(Udep_))

# k = m - matrix_rank(E)

# factorization = scipy.linalg.lu(E.transpose())
# P, L, U = factorization

# print('E shape/rank:', E.shape, matrix_rank(E))
# print('P shape/rank:', P.shape, matrix_rank(P))
# print('L shape/rank:', L.shape, matrix_rank(L))
# print('U shape/rank:', U.shape, matrix_rank(U))
# print('k:', k)

# L1 = L[:m - k]
# L2 = L[m - k:]

# print('L1 shape', L1.shape)
# print('L2 shape', L2.shape)

# LTinv = pinv(L1).transpose()
# UTinv = pinv(U).transpose()

# d = zeros((m, 1))
# d[:m - k] = LTinv.dot(UTinv).dot(ph)

# B = zeros((m, k))
# Ik = identity(k)
# L1L2 = - LTinv.dot(L2.transpose())
# print(L1L2.shape)
# print(P.shape)

# B[:m - k] = L1L2
# B[m - k:] = Ik

# B = P.dot(B)

# # print(M.ind)

# qid = ones((k, 1))

# q = B.dot(qid) + d

# print(q)

# print(P.dot(qid))

# print('max(q), min(q)', max(q), min(q))

# print('max(ph), min(ph)', max(ph), min(ph))

# Eq = E.dot(q)

# print('max(Eq), min(Eq)', max(Eq), min(Eq))

# print('B:', B.shape)
# print('d:', d.shape)

# q = M.q
# print('q shape:', q.shape)

# qid = q[M.ind]

# print(qid)

# q2 = B.dot(qid) + d

# print('q2 shape:', q2.shape)

# diffq = abs(q - q2)**2

# print(sum(diffq))

# q = B.dot(qid) + d

# print('ind:', len(ind))

