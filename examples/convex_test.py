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
from compas_tno.utilities import check_constraints
from compas_tno.utilities import circular_heights
from compas_tno.utilities.symmetry import not_sym_load
from compas_tno.utilities import replicate_contraints

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

def save_matlab(form, file, find_inds=True, heights=False, lb_ub=True):

    args = initialize_problem(form, find_inds=find_inds)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args

    if lb_ub:
        zt = []
        zlb = []
        zub = []

        for key in form.vertices():
            if form.get_vertex_attribute(key, 'is_fixed') == False:
                zt.append(form.get_vertex_attribute(key,'target'))
                zlb.append(form.get_vertex_attribute(key,'lb'))
                zub.append(form.get_vertex_attribute(key,'ub'))
    else:
        zt = zlb = zub = 1.0

savemat(file,{
    'm': int(form.number_of_edges()),
    'n': int(form.number_of_vertices()),
    'ni': len(free),
    'C': C,
    'Ci': Ci,
    'Cf': Cf,
    'Cit': Cit,
    'x': x,
    'y': y,
    'xb': x[fixed],
    'yb': y[fixed],
    'pz': pz[free],
    'free': free,
    'fixed': fixed,
    'p': p,
    'E': E,
    'Edinv': Edinv,
    'Ei': Ei,
    'k': k,
    'ind': ind,
    'dep': dep,
    'zlb': zlb,
    'zub': zub,
    'zt': zt
})
    pass


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    j = 2
    i = 6

    i_s = []
    n_s = []
    lp_s = []
    calc = []

    # ['B1','B2','B3','C1','C2','C3','C4','C5','D1','D2','D3','D4','D5']

    for i in ['B1']:
        print(i)
        j = 2
        i = 6

        # Load

        # file = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/01.json'
        file_constraints = '/Users/mricardo/compas_dev/me/loadpath/Fix/nosym/0'+str(j)+'_0'+str(i)+'_t_60.json'
        file = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_complete_nosym.json'
        file = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_complete_paper.json'

        file_save = '/Users/mricardo/compas_dev/me/loadpath/Fix/nosym/0'+str(j)+'_0'+str(i)+'_t_60cm_p_120.json'
        # file = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0'+str(j)+'_0'+str(i)+'_complete_paper.json'
        # file_arch = '/Users/mricardo/compas_dev/me/minmax/2D_arch/01_lp.json'

        # file = '/Users/mricardo/compas_dev/me/convex/4bars/diagram.json'
        # file_matlab_solve = '/Users/mricardo/Documents/MATLAB/optimisation/topology/Corner_0'+str(j)+'_0'+str(i)+'_solve.mat'
        # file_matlab = '/Users/mricardo/Documents/MATLAB/optimisation/Corner_0'+str(j)+'_0'+str(i)+'.mat'
        file_matlab = '/Users/mricardo/Documents/MATLAB/optimisation/discretize/nosym_02_06_p_140.mat'
        # file_pic = '/Users/mricardo/Documents/MATLAB/optimisation/topology/Corner'+str(i)+'_pic.pdf'

        form = FormDiagram.from_json(file)
        # form = _form(form)
        # form = create_sym(form)

        # Modify Form

        # form = fix_boundaries_complete(form)
        # form = fix_mid_complete(form)
        # form = set_cross_vault_heights(form, ub_lb=True, thk=2.0, set_heights=False)
        # form = circular_heights(form, thk = 0.20)
        form = replicate_contraints(file,file_constraints)
        form = not_sym_load(form, magnitude = 1.40)
        check_constraints(form, show=True)
        overview_forces(form)
        plot_form(form,show_q=False, max_width=2.0, simple=True).show()
        # form.to_json(file_fix)

        # Save Matlab

        save_matlab(form, file_matlab, find_inds=True)
        plot_form(form).show()

        # Load Matlab

        form = load_matlab(form,file_matlab)
        if form.attributes['loadpath'] < 1000.0:
            overview_forces(form)
            i_s.append(i)
            n_s.append(form.number_of_edges())
            lp_s.append(form.attributes['loadpath'])
            if form.attributes['loadpath'] < 1000.0:
                calc.append(i)
                plot_form(form,show_q=False, max_width=2.0, simple=True).show()
                form.to_json(file_save)
                check_constraints(form, show=True)

    print(i_s)
    print(n_s)
    print(lp_s)
    print(calc)

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()
