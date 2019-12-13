from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
# from compas_tno.algorithms.equilibrium import reactions
from compas_tno.plotters.plotters import plot_form
from compas_tno.plotters.plotters import plot_force
from compas_tno.algorithms import initialize_problem
from compas_tno.algorithms import update_tna
from compas_tno.algorithms import z_from_form
from compas_tno.diagrams.form import overview_forces
from compas.utilities import geometric_key

from compas_tno.utilities import fix_boundaries_sym
from compas_tno.utilities import fix_boundaries_complete
from compas_tno.utilities import set_cross_vault_heights
from compas_tno.utilities import set_pavillion_vault_heights
from compas_tno.utilities import set_oct_vault_heights

from compas_tno.utilities import fix_mid_sym
from compas_tno.utilities import fix_mid_complete

from compas.geometry import is_point_on_segment
from compas.geometry import intersection_segment_segment
from compas.geometry import is_intersection_line_line
from compas.geometry import intersection_line_line
from numpy import shape
from numpy import abs
from numpy import argmin
from numpy import array
from numpy import float64
from numpy import dot
from numpy import hstack
from numpy import isnan
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum
from numpy import vstack
from numpy import zeros
from numpy import ones
from numpy import append
# from numpy import todense
from numpy.linalg import pinv
from numpy.linalg import matrix_rank
from numpy.random import rand
from numpy.random import randint

from compas_tno.diagrams.form import _form
from compas_tno.algorithms import min_loadpath
from compas_tno.algorithms import min_thrust

from compas.numerical import normrow
from compas_viewers.meshviewer import MeshViewer
# from compas.numerical import norm

from scipy.linalg import svd
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from compas.numerical import equilibrium_matrix

from compas_tno.utilities.constraints import set_cross_vault_heights
from compas_tno.utilities.constraints import set_pavillion_vault_heights
from compas_tno.utilities.constraints import circular_heights
from compas_tno.utilities.symmetry import create_sym

import matplotlib
import matplotlib.pyplot as plt

from scipy.io import savemat
from scipy.io import loadmat

def load_matlab(form,file):

    sol = loadmat(file)
    q = sol['q']
    cvx_status = sol['cvx_status']
    print(cvx_status)

    if cvx_status is not (['Infeasible'] or ['Failed']):
        uv_i = form.uv_index()
        for u,v in form.edges():
            i = uv_i[(u,v)]
            [qi]= q[i]
            # print(qi)
            form.set_edge_attribute((u,v), 'q', qi)

        form = z_from_form(form)
        return form
    else:
        print('Failed')

    return None

def save_matlab(form, file, find_inds=True, heights=False):

    args = initialize_problem(form, find_inds=find_inds)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    w = C.dot(z)
    Ez = Cit.dot(diags(w))

    savemat(file,{
            'm': int(form.number_of_edges()),
            'n': int(form.number_of_vertices()),
            'C': C,
            'Ci': Ci,
            'Cf': Cf,
            'Cit': Cit,
            'xt': x.transpose(),
            'yt': y.transpose(),
            'xf': x[fixed],
            'yf': y[fixed],
            'pz': pz[free],
            'p': p,
            'E': E,
            'Edinv': Edinv,
            'Ei': Ei,
            'k': k,
            'ind': ind,
            'dep': dep,
            'Ez': Ez
        })
    pass


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Load

    # import sdpt3glue

    i = 4
    j = 2
    # file = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0'+str(j)+'_0'+str(i)+'_complete_paper.json'
    # file_matlab = '/Users/mricardo/Documents/MATLAB/optimisation/discretize/Corner_0'+str(j)+'_0'+str(i)+'.mat'
    file_matlab = '/Users/mricardo/Documents/MATLAB/optimisation/2Darch_solve.mat'
    file = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/01.json'
    form = FormDiagram.from_json(file)
    plot_form(form,heights=True, show_q=False).show()
    # form = circular_heights(form)

    # form = _form(form)
    # form = create_sym(form)

    # Modify Form

    # form = fix_boundaries_complete(form)
    # form = fix_mid_complete(form)
    # plot_form(form,show_q=False, max_width=2.0).show()
    # form.to_json(file_fix)

    # Save Matlab

    # save_matlab(form, file_matlab, find_inds=False)

    # Load Matlab

    # form = load_matlab(form, file_matlab)

    # Run CVOPT

    # args = initialize_problem(form, find_inds=False)
    # q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args

    # form = min_loadpath(form, args, printout=True)

    # zmin = zeros(len(free)).reshape(len(free),1)
    # zmax = 2.0*ones(len(free)).reshape(len(free),1)
    # print(zmin.shape,zmax.shape)

    # form = min_thrust(form, args, zmin, zmax, printout=True)

    # solver = form.attributes['solve']
    # print('Solver utilised: ', solver.solver_name)
    # print('Solve time: {0}'.format(solver.solve_time))
    # print('Set-up time: {0}'.format(solver.setup_time))

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()
