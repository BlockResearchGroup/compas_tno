from compas_tna.diagrams import FormDiagram

from compas_tno.utilities.constraints import circular_heights
from compas_tno.utilities.constraints import circular_joints
from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_arch
from compas_tno.algorithms import optimise_general

from compas_tno.plotters.plotters import plot_form_xz
from compas_tno.plotters.plotters import plot_form_joints
from numpy import array

from compas.geometry import intersection_segment_segment_xy
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import intersection_line_line_xy

from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import divide
from numpy import array
from numpy.linalg import matrix_rank
from numpy import transpose
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import scale_vector_xy
from compas.geometry import distance_point_point_xy
from compas_tno.algorithms import zlq_from_qid
from compas_tno.algorithms import q_from_qid
from scipy.sparse import diags

from compas_tno.algorithms.problems import initialise_problem

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file_save = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/arch_joint.json'
    lambda_hor = 0.0
    exitflag = 0
    blocks = 18
    form = create_arch(total_nodes = blocks, lambda_hor = lambda_hor, D = 2.00)
    # form.to_json(file_save)
    print_opt = True
    thk = 0.25
    form = circular_joints(form, thk = thk, blocks = blocks)
    # plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()

    # Initial parameters

    translation = form.attributes['tmax']
    bounds_width = 5.0
    use_bounds = False
    qmax = 10000
    indset = None

    # Pre-Otimisation

    exitflag = 1
    count = 0

    while count < 10 and exitflag is not 0:
        # form = _form(form, keep_q=True)
        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver='slsqp',
                                            printout=print_opt,
                                            qmin = 0.1,
                                            find_inds=True,
                                            tol=0.01,
                                            translation = translation,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='constr_lp',
                                            indset=indset,
                                            bmax = True,
                                            summary=print_opt )

        # Check compression and Save

        q = [attr['q'] for u, v, attr in form.edges(True)]
        qmin  = min(array(q))
        count += 1
        if qmin > -0.1 and exitflag == 0:
            print('Optimisation completed - Trial:',count, 'Lambda Px', lambda_hor)
            # plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()

    # Optimisation

    exitflag = 1
    count = 0

    while count < 10 and exitflag is not 0:
        # form = _form(form, keep_q=True)
        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver='slsqp',
                                            printout=print_opt,
                                            qmin = 0.01,
                                            find_inds=True,
                                            tol=0.01,
                                            translation = translation,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='min',
                                            indset=indset,
                                            bmax = True,
                                            summary=print_opt)

        # Check compression and Save

        q = [attr['q'] for u, v, attr in form.edges(True)]
        qmin  = min(array(q))
        count += 1
        if qmin > -0.1 and exitflag == 0:
            print('Optimisation completed - Trial:',count, 'Lambda Px', lambda_hor)
            # plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()
            # form.to_json(file_save)


    # Optimisation

    exitflag = 1
    count = 0
    qmax = qopt[0]*1.1
    qmin = qopt[0]*0.9

    while count < 100:
        # form = _form(form, keep_q=True)
        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver='ga',
                                            printout=print_opt,
                                            qmin = qmin,
                                            find_inds=True,
                                            tol=0.01,
                                            translation = translation,
                                            tension=False,
                                            use_bounds = use_bounds,
                                            bounds_width = bounds_width,
                                            objective='min_joints',
                                            indset=indset,
                                            bmax = True,
                                            summary=print_opt)

        # Check compression and Save

        q = [attr['q'] for u, v, attr in form.edges(True)]
        qmin  = min(array(q))
        count += 1
        if qmin > -0.1 and exitflag == 0:
            print('Optimisation completed - Trial:',count, 'Lambda Px', lambda_hor)
            plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()
            # form.to_json(file_save)

    # for i, elements in joints.items():
    #     p1, p2, edges = elements
    #     for k1, k2 in edges:
    #         if k1 in fixed:
    #             joints[i][2].add((-k1,k1))
    #         if k2 in fixed:
    #             joints[i][2].add((-k2,k2))
    # print(joints)

    # i_k = form.index_key()
    # result = []
    # for i, elements in joints.items():
    #     p1, p2, edges = elements
    #     # print('Joint: ',i)
    #     joint_int = - 50.0
    #     for k1, k2 in edges:
    #         # print('Edge',k1,k2)
    #         if -k1 == k2:
    #             e1 = xR[-k1], zR[-k1]
    #         else:
    #             e1 = form.vertex_coordinates(i_k[k1])[0], form.vertex_coordinates(i_k[k1])[2]
    #         e2 = form.vertex_coordinates(i_k[k2])[0], form.vertex_coordinates(i_k[k2])[2]
    #         edge_2D = [e1,e2]
    #         joint_2D = [[p1[0],p1[2]],[p2[0],p2[2]]]
    #         # print(is_intersection_segment_segment_xy(joint_2D,edge_2D))
    #         pt = intersection_line_line_xy(edge_2D,joint_2D)
    #         # print(pt)
    #         if is_intersection_segment_segment_xy(joint_2D,edge_2D) == True:
    #             joint_int = 0.0
    #             break
    #         else:
    #             virtual_pt = intersection_line_segment_xy(joint_2D,edge_2D)
    #             # print(virtual_pt)
    #             if virtual_pt:
    #                 offset = min(distance_point_point_xy([virtual_pt[0],virtual_pt[2]],joint_2D[0]),distance_point_point_xy([virtual_pt[0],virtual_pt[2]],joint_2D[1]))
    #                 # print(offset)
    #                 joint_int = -offset
    #                 break
    #     # print('Joint Int',joint_int)
    #     result.append(joint_int)

    # print(result)

    plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()
