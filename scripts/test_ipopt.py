from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form

from compas_tno.utilities.constraints import set_cross_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.algorithms.general_solver import set_b_constraint
from compas_tno.algorithms.general_solver import set_joints_constraint
from compas_tno.algorithms.general_solver import set_cracks_constraint

from compas_tno.plotters.plotters import plot_form

from scipy.optimize import fmin_slsqp
from compas.numerical import devo_numpy
from compas.numerical import ga

from compas_tno.algorithms.problems import initialise_problem

from compas_tno.algorithms.objectives import f_min_loadpath
from compas_tno.algorithms.objectives import f_min_loadpath_pen
from compas_tno.algorithms.objectives import f_min_thrust
from compas_tno.algorithms.objectives import f_min_thrust_pen
from compas_tno.algorithms.objectives import f_max_thrust
from compas_tno.algorithms.objectives import f_target
from compas_tno.algorithms.objectives import f_constant

from compas_tno.algorithms.constraints import f_compression
from compas_tno.algorithms.constraints import f_ub_lb
from compas_tno.algorithms.constraints import f_joints
from compas_tno.algorithms.constraints import f_cracks

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import zlq_from_qid
from compas.utilities import geometric_key

import pyOpt

from numpy.random import rand
from numpy import append
from numpy import array

from compas_tna.diagrams import FormDiagram

import math

# import ipopt

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

type_fd = 'cross_fd'
objective = 'min'
thck = 0.50

x_span = 5.0
y_span = 10.0
example = 'rectangular/5x10/'

if type_fd == 'cross_fd':
    divisions = 20
    form = create_cross_form(xy_span=[[0.0, x_span], [0.0, y_span]], division=divisions)

if type_fd == 'fan_fd':
    divisions = 16
    form = create_fan_form(xy_span=[[0.0, x_span], [0.0, y_span]], division=divisions)

PATH = '/Users/mricardo/compas_dev/me/minmax/cross/' + example + type_fd + '/' + type_fd + '_discr_' + str(divisions)
file_initial = PATH + '_lp.json'
# file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'

form = FormDiagram.from_json(file_initial)
indset = form.attributes['indset']

translation = True
bounds_width = 5.0
use_bounds = False
qmax = 100
qmin = -1e-6
print_opt = True

form = set_cross_vault_heights(form, xy_span=[[0.0, x_span], [0.0, y_span]], thk=thck, b=5.0, set_heights=False, ub_lb=True, update_loads=True)

k_i = form.key_index()
i_k = form.index_key()
i_uv = form.index_uv()
printout = True

find_inds = True
tol = 0.0001
args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds, tol=tol)
q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

bmax = True
cracks = False
b = set_b_constraint(form, bmax, printout)
joints = set_joints_constraint(form, printout)
cracks_lb, cracks_ub = set_cracks_constraint(form, cracks, printout)

args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty)

fobj, fconstr = f_min_thrust, f_ub_lb

if translation:
    x0 = q[ind]
    zb_bounds = [[form.get_vertex_attribute(i_k[i], 'lb'), form.get_vertex_attribute(i_k[i], 'ub')] for i in fixed]
    bounds = [[qmin, qmax]] * k + zb_bounds
    x0 = append(x0, z[fixed]).reshape(-1, 1)
else:
    x0 = q[ind]
    bounds = [[qmin, qmax]] * k

print('Total of Independents:', len(ind))
print('Number of Variables:', len(x0))
f0 = fobj(x0, *args)
g0 = fconstr(x0, *args)
print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))

from ipopt import minimize_ipopt

minimize_ipopt(fobj, x0, args = args, constraints = [fconstr])

