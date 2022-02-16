import json
from numpy import array
from numpy import zeros

import cvxpy as cp
from cvxpy import diag
from cvxpy import matrix_frac
from cvxpy import Minimize
from cvxpy import Problem

path = '/Users/mricardo/compas_dev/compas_tno/data/data.json'

with open(path, 'r') as json_file:
    data = json.load(json_file)

q = array(data['q'])
E = array(data['E'])
C = array(data['C'])
Ci = array(data['Ci'])
Cit = array(data['Cit'])
Cb = array(data['Cb'])
pz = array(data['pz'])
qmin = array(data['qmin'])
qmax = array(data['qmax'])
free = data['free']
fixed = data['fixed']
x = array(data['x'])
y = array(data['y'])
m = data['m']

import scipy

factorization = scipy.linalg.lu(E)
print(factorization)
P, L, U = factorization

print('E:', E.shape)
print('P:', E.shape)
print('L:', E.shape)
print('U:', E.shape)
# print('ind:', len(ind))


# # print(qmin)  # -100000
# # print(qmax)  # 1e-8

# q = cp.Variable(m)
# scale = 1.0

# # fobj = matrix_frac(pz[free], Cit@cp.diag(q)@Ci) + x.T@C.T@diag(q)@Cb@x[fixed] + y.T@C.T@diag(q)@Cb@y[fixed]
# fobj = matrix_frac(pz[free], - scale*Cit@cp.diag(q)@Ci) - scale*x.T@C.T@diag(q)@Cb@x[fixed] - scale*y.T@C.T@diag(q)@Cb@y[fixed]
# objective = Minimize(fobj)

# horz = E@q == 0
# # pos = q >= zeros(m)
# # maxq = q <= -1/100 * qmin.flatten()
# # pos = q >= 2000.0
# # maxq = q <= 0
# pos = q >= qmin.flatten()
# maxq = q <= qmax.flatten()

# constraints = [horz, pos, maxq]

# prob = Problem(objective, constraints)
# prob.solve(verbose=True, solver='MOSEK')

# q = q.value

# print('max / min pz:', max(pz[free]), min(pz[free]))
# print('max / min q:', max(q), min(q))
# print('max / min qmin:', max(qmin), min(qmin))
# print('max / min qmax:', max(qmax), min(qmax))

# # fobj_1 = matrix_frac(pz[free], - Cit@cp.diag(q)@Ci)
# # fobj_2 = - x.T@C.T@diag(q)@Cb@x[fixed]
# # fobj_3 = - y.T@C.T@diag(q)@Cb@y[fixed]

