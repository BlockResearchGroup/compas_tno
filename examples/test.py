from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
# from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.plotters.plotters import plot_form
from compas_thrust.plotters.plotters import plot_force
from compas_thrust.algorithms import initialize_problem
from compas_thrust.algorithms import update_tna
from compas_thrust.algorithms import z_from_form
from compas_thrust.diagrams.form import overview_forces
from compas.utilities import geometric_key

from compas_thrust.utilities import fix_boundaries_sym
from compas_thrust.utilities import fix_boundaries_complete
from compas_thrust.utilities import set_cross_vault_heights
from compas_thrust.utilities import set_pavillion_vault_heights
from compas_thrust.utilities import set_oct_vault_heights

from compas_thrust.utilities import fix_mid_sym
from compas_thrust.utilities import fix_mid_complete

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
from numpy.linalg import pinv
from numpy.linalg import matrix_rank
from numpy.random import rand
from numpy.random import randint

from compas.numerical import normrow
from compas_viewers.meshviewer import MeshViewer
# from compas.numerical import norm

from scipy.linalg import svd
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

import matplotlib
import matplotlib.pyplot as plt


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # i=3
    # filesym = '/Users/mricardo/compas_dev/me/minmax/radial/02_0'+ str(i) +'_calc.json'
    # form = FormDiagram.from_json(filesym)
    # # ind = find_independents(E)
    # uv_i = form.uv_index()
    # form.update_default_edge_attributes({ 'is_symmetry': False })
    
    # args = initialize_problem(form)
    # q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    # print(sym)
    # print(ind)
    # for u,v in form.edges():
    #     if form.get_edge_attribute((u,v),'is_ind') == True:
    #         print(uv_i[(u,v)])
    #         print(form.get_edge_attribute((u,v), 'q'))
    #         print(form.get_edge_attribute((u,v),'is_symmetry'))

    # print('next to try')
    # print(uv_i[(10,0)])

    # plot_form(form,show_edgeuv=True, max_width=2.0).show()

    i_s = []
    n_s = []
    lp_s = []

    for i in range(5,6):

        print('\n\n------------------ Form ',str(i),'\n')
        j = 1

        file_complete = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+ str(j) +'_0'+ str(i) +'_complete.json'
        form = FormDiagram.from_json(file_complete)
        overview_forces(form)
        form = set_cross_vault_heights(form, set_heights=True)
        plot_form(form, heights=True, show_q = False).show()
        # plot_form(form).show()
        # form = z_from_form(form)
        # plot_form(form).show()
        # form.to_json(file_complete)
        # i_s.append(i)
        # n_s.append(form.number_of_edges())
        # lp_s.append(form.attributes['loadpath'])
        # plot_form(form,show_q=False, simple=True, max_width=5.0).show()
        # force = ForceDiagram.from_formdiagram(form)
        # form, force = update_tna(form, delete_face=False)
        # plot_force(force, form, color_inds=False).show()

        # form = not_sym_load(form, 0.0, 5.0, 25)

        viewer = MeshViewer()
        viewer.mesh = form
        viewer.show()



    # fig, ax = plt.subplots()
    # ax.plot(n_s, lp_s)
    # ax.set(xlabel='# Edges F.D.', ylabel='LP',
    #     title='LP for different discretization - Orthogonal')
    # ax.grid()
    # plt.show()

    # print(lp_s)


    # for i in range(1,8):
    #     print('Form: ',str(i))
    #     filesym = '/Users/mricardo/compas_dev/me/minmax/radial/01_0'+ str(i) +'_calc.json'
    #     form_ = FormDiagram.from_json(filesym)
    #     args = initialize_problem(form_)
    #     plot_form(form_).show()
    #     file = '/Users/mricardo/compas_dev/me/minmax/radial/01_0'+ str(i) +'_complete_min.json'
    #     form = FormDiagram.from_json(file)

    #     inds_midpt = []
    #     for u,v in form_.edges():
    #         if form_.get_edge_attribute((u,v),'is_ind') == True:
    #             inds_midpt.append(geometric_key(form_.edge_midpoint(u,v)[:2]+[0]))
    #     for u,v in form.edges():
    #         form.set_edge_attribute((u,v), 'is_ind', value = False)
    #         if geometric_key(form.edge_midpoint(u,v)) in inds_midpt:
    #             print('update ',u,v)
    #             form.set_edge_attribute((u,v), 'is_ind', value = True)
        
    #     for u,v in form.edges():
    #         if geometric_key(form.edge_midpoint(u,v)[:2]+[0]) in inds_midpt:
    #             print('update ',u,v)
    #             form.set_edge_attribute((u,v), 'is_ind', value = True)
                
    #     plot_form(form,radius = 0.05, show_q= False, max_width=2).show()
    #     form, force = update_tna(form)
    #     plot_force(force, form).show()

        # form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
        # form = z_from_form(form)
        # plot_form(form,radius = 0.05, show_q= False, fix_width=True, max_width=5).show()
        # overview_forces(form)
        # file_out = '/Users/mricardo/compas_dev/me/minmax/radial/02_0'+ str(i) +'_complete_fdm.json'
        # form.to_json(file_out)
        # form, force = update_tna(form)
        # plot_force(force, form).show()
        # plot_form(form, max_width = 5, fix_width=False, show_q= False).show()
        # reactions(form, plot=True)
    #     # form.to_json(file)

    #     # lines = [
    #     #         [   [0.0,0.0,0.0],[0.0,2.0,0.0]     ],
    #     #         [   [0.0,0.0,0.0],[0.0,-2.0,0.0]    ],
    #     #         [   [0.0,0.0,0.0],[2.0,0.0,0.0]     ],
    #     #         [   [0.0,0.0,0.0],[-2.0,0.0,0.0]    ]
    #     #         ]
    #     # form = FormDiagram.from_lines(lines, delete_boundary_face=False)
    #     # form.update_default_vertex_attributes({'is_roller': False})
    #     # form.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    #     # print(form.number_of_edges())
    #     # for key in form.vertices():
    #     #     if form.vertex_coordinates(key) == [0,0,0]:
    #     #         pass
    #     #     else:
    #     #         form.set_vertex_attribute(key, 'is_fixed', True)
    #     #         form.set_vertex_attribute(key, 'is_anchored', True)
    #     # form.plot()
    #     # plot_form(form, fix_width= 0.1, simple=True).show()
    #     # force = ForceDiagram.from_formdiagram(form)
    #     # plot_force(force,form)

    # file = '/Users/mricardo/compas_dev/me/minmax/radial/02_02_complete.json'
    # form = FormDiagram.from_json(file)
    # print(form.number_of_edges())
    # form.plot()
    # plot_form(form, fix_width=True, max_width=5, show_q=False, simple=True).show()
    # plot_form(form).show()
    # q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = initialize_problem(form)
    # print(E.shape)
    # CfQ = Cf.transpose().dot(diags(q.flatten()))
    # rx = (CfQ.dot(U[:,newaxis]))
    # ry = (CfQ.dot(V[:,newaxis]))
    # uv_i = form.uv_index()
    # print(rx)
    # print(ry)
    # for u, v in form.edges():
    #     form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)
    # print(norm(rx))
    # print(norm(ry))
    # print(ind)
    # C = C.todense()

    # plot_form(form, fix_width=True, max_width=5, show_q=False, simple=False).show()
    # form, force = update_tna(form, delete_face = True)
    # force = ForceDiagram.from_formdiagram(form)
    # force.plot()

    # for i in range(C.shape[0]):
    #     print(C[i,])