from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.grad_based import optimise_tna
from compas_tno.algorithms.scale import evaluate_scale
from compas_tno.algorithms.scale import lagrangian_scale
from compas_tno.algorithms.scale import scale_form
from compas_tno.algorithms.equilibrium import z_from_form

from compas_tno.algorithms import planes_trimesh
from compas_tno.algorithms import local_matrix
from compas_tno.algorithms import local_matrix_external
from compas_tno.algorithms import assembly_Cf
from compas_tno.algorithms import A_heights
from compas_tno.algorithms import A_stress
from compas_tno.algorithms import hessian
from compas_tno.algorithms import simple_nurbs

from compas_tno.utilities.constraints import check_constraints
from compas_tno.utilities.symmetry import replicate
from compas_tno.utilities import paraboloid
from compas_tno.utilities import dome

from compas_tno.diagrams.form import energy
from compas_tno.diagrams.form import loadpath
from compas_tno.diagrams.form import adapt_tna
from compas_tno.diagrams.form import evaluate_a
from compas.geometry import Plane


from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_dual

# from scipy import array
# from scipy import tensordot
# from scipy.optimize import nnls
# from scipy import tensorproduct

from sympy import Array
from sympy import tensorproduct

from compas.numerical import normrow
from compas.numerical import normalizerow
from compas.geometry import normalize_vector
from compas.datastructures import mesh_face_matrix
from compas.datastructures import mesh_quads_to_triangles
from compas.geometry import is_ccw_xy

from compas.datastructures import Mesh
from compas_plotters import MeshPlotter

# from compas.geometry import Mesh
from compas_viewers.meshviewer import MeshViewer

from copy import deepcopy
from numpy import array
from numpy import argmin
from numpy import cross
from numpy import dot
from numpy import zeros
from numpy import int8
from numpy.linalg import lstsq
from numpy.linalg import matrix_rank
from numpy.linalg import pinv
from numpy.linalg import det
from numpy import delete
from numpy import vstack
from numpy import hstack

from compas.numerical.linalg import spsolve
from compas.numerical.linalg import spsolve_with_known
from compas.numerical.linalg import solve_with_known

import math
from numpy.linalg import inv

def view_form(form, plot_z = None):

    if plot_z is not None:
        for key in form.vertices():
            zi = float(plot_z[k_i[key]])
            form.vertex_attribute(key, 'z', zi)

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()

    return

def update_heights(form, f):

    for key in form.vertices():
        form.vertex_attribute(key, 'z', value=float(f[k_i[key]]))

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/me/lsm/form/6node_3.json'
    # file = '/Users/mricardo/compas_dev/me/loadpath/Fix/A_lp.json'
    # file = '/Users/mricardo/compas_dev/me/discretize/01_04_complete.json'

    # surf = simple_nurbs(5,5,10, plot=True)

    form = FormDiagram.from_json(file)
    form.delete_face(0)
    q = []
    Rm_ = []

    plot_form(form).show()
    form = z_from_form(form)

    # mesh_quads_to_triangles(form)
    for u,v in form.edges():
        qi = form.edge_attribute((u,v),'q')
        # if qi == 1.0 or qi == None:
        #     form.edge_attribute((u,v),'q',value=0.0)
        Rm_.append(form.edge_attribute((u,v),'q'))

    plot_form(form, fix_width=True, max_width=2.0, show_q=True).show()
    print('Form Has {0} edges and {1} nodes'.format(form.number_of_edges(),form.number_of_vertices()))

    k_i = form.key_index()
    i_k = form.index_key()
    uv_i = form.uv_index()
    i_uv = form.index_uv()

    # Apply shape function f

    V = []
    phi = []
    index_bound = []
    nodes_bound = []
    E_bound = form.edges_on_boundary()
    N_bound = form.vertices_on_boundary()
    for pair in E_bound:
        u,v = pair
        try:
            index_bound.append(uv_i[(u,v)])
        except BaseException:
            index_bound.append(uv_i[(v,u)])
    zs = []

    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if xi == 0.0 and yi == 0.0:
            key_center = k_i[key]
        zi = paraboloid(xi,yi)
        # zi = dome(xi,yi,10.1)
        phii = zi
        # zi = surf.evaluate_single(((xi+5)/10,(yi+5)/10))[2]
        q.append(form.vertex_attribute(key,'pz'))
        if key in N_bound:
            nodes_bound.append(k_i[key])
        # q.append(0.0)
        # form.vertex_attribute(key, 'z', zi)
        # form.vertex_attribute(key, 'phi', phii)
        # phi.append(phii)
        zs.append(form.vertex_attribute(key,'z'))

    # view_form(form)
    print('zs: {0:.3f} : {1:.3f}'.format(float(min(zs)), float(max(zs))))
    plot_form(form,show_q=False,heights = True).show()

    # view_form(form)

    # ---------------- 1. Assembly Matrix

    Cf = assembly_Cf(form)

    # --------------- 2. Get Stress Function From zi's

    A = A_heights(form)

    # ---------------- 2.1. Request Rm > 0 ( Convex )

    # print('Matrix A Shape:',A.shape)
    # Rm, res = nnls(A, array(q))
    # print('Rm: {0:.3f} : {1:.3f}'.format(float(min(Rm)), float(max(Rm))))
    # print('Residual Rm:', res)
    # print('Rm shape:', Rm.shape)

    # sol = lstsq(Cf,Rm,rcond=None)
    # phi = sol[0]
    # res = sol[1]

    # print('phi: {0:.3f} : {1:.3f}'.format(float(min(phi)), float(max(phi))))
    # print('Residual Phi', res)
    # print('Phi shape:', phi.shape)

    # ---------------- 2.1. Not request Rm > 0 ( Convex )

    # K = dot(A, Cf)
    # print('Matrix K Shape:',K.shape)
    # print('Matrix K Det:',det(K))
    # print('Matrix K Rank:',matrix_rank(K))
    # # print(K)
    # sol = lstsq(K,array(q).reshape(len(q),1),rcond=None)
    # phi = sol[0]
    # res = sol[1]
    # print('phi: {0:.3f} : {1:.3f}'.format(float(min(phi)), float(max(phi))))
    # print('residual phi (Not Concave Constraint)', res)
    # print(phi.shape)

    Rm = array(Rm_).reshape(len(Rm_),1)
    # Rm = dot(Cf,phi)
    print(Rm)
    print(Rm.shape)
    print('Rm: {0:.3f} : {1:.3f}'.format(float(min(Rm)), float(max(Rm))))

    # ------------- 3. Get Z's from Stress Function Phi

    form.update_default_edge_attributes({'Rm' : 0.0}) # Change to RM attribute
    for u,v in form.edges():
        form.edge_attribute((u,v), 'Rm', value= float(Rm[uv_i[(u,v)]]))

    plot_form(form, thick='Rm', show_q= False, heights=True, radius=0.3).show()

    A = A_stress(form,Rm)
    # Check if it is the same as [ Ct (Rm/lij) C ]

    # Known Var

    # A = delete(A, key_center, axis=0)
    # A = delete(A, key_center, axis=1)
    # q = delete(q, key_center)

    print('Matrix A Shape:',A.shape)
    print('Matrix A Det:',det(A))
    print('Matrix A Rank:',matrix_rank(A))
    # print(K)
    # sol = lstsq(A,array(q).reshape(len(q),1),rcond=None)
    sol = solve_with_known(A,array(q).reshape(len(q),1),zeros((len(q),1)),nodes_bound)
    # sol = spsolve(A,array(q).reshape(len(q),1))#,rcond=None)
    f = sol
    res = 0
    # f = sol[0]
    # res = sol[1]
    # f = f-1270.0
    # f[f<0] = 0.0
    # f = f/2
    # f = f -654800
    # f[f<0] = 0.0
    print('f: {0:.3f} : {1:.3f}'.format(float(min(f)), float(max(f))))
    print('Residual f', res)
    # print(f)
    # f = vstack([f[:key_center],[0.0],f[key_center:]])



    # --------------- 2. Get Stress Function From zi's

    update_heights(form,f)
    A = A_heights(form)
    K = dot(A, Cf)
    print('Matrix K Shape:',K.shape)
    print('Matrix K Det:',det(K))
    print('Matrix K Rank:',matrix_rank(K))
    # print(K)
    sol = lstsq(K,array(q).reshape(len(q),1),rcond=None)
    phi = sol[0]
    res = sol[1]
    print('phi: {0:.3f} : {1:.3f}'.format(float(min(phi)), float(max(phi))))
    print('residual phi (Not Concave Constraint)', res)
    print(phi.shape)

    # -------------- 4. Plot Results

    plot_form(form,show_q=False, heights = True).show()
    view_form(form, f)
    # view_form(form, phi)

    plot_dual(form).show()
