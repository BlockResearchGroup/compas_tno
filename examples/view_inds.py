
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single
from compas_thrust.algorithms.ind_based import initialize_problem

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check
from compas_thrust.algorithms.equilibrium import zlq_from_qid
from compas_thrust.algorithms.equilibrium import update_qid

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.diagrams.form import oveview_forces
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import _form
from compas_thrust.plotters.plotters import plot_form

from compas_thrust.plotters.plotters import plot_form
# from compas_viewers.meshviewer import MeshViewer
from compas_thrust.viewers.meshviewer import MeshViewer

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":


    file = '/Users/mricardo/compas_dev/me/minmax/2D_arch/02_lp.json'
    # file_complete = '/Users/mricardo/compas_dev/me/minmax/radial/02_06_complete.json'

    form = FormDiagram.from_json(file)
    form = update_qid(form, 100.0)

    # Create a Function to do that

    # args = initialize_problem(form)

    # q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    # k_i = form.key_index()
    # i_uv = form.index_uv()
    # ind = args[1]
    
    # print(ind)
    
    # q0 = []
    # for i in ind:
    #     key = i_uv[i]
    #     q0.append(form.get_edge_attribute(key, 'q'))
    
    # print(q0)

    # # Modify via Sliding
    # q0[0] += 100.0

    # # Recalculate and Update
    # q[ind, 0] = q0
    # q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
    # z[free] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free])

    # for key, attr in form.vertices(True):
    #     index = k_i[key]
    #     attr['z']  = z[index]

    # View

    viewer = MeshViewer()
    # viewer.mesh = form
    # viewer.setup()
    viewer.show()
