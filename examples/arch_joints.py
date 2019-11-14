
from compas_tna.diagrams import FormDiagram
from compas_thrust.algorithms import optimise_general

from compas_thrust.utilities.constraints import circular_heights
from compas_thrust.utilities.constraints import circular_joints
from compas_thrust.diagrams.form import overview_forces
from compas_thrust.diagrams.form import create_arch

from compas_thrust.plotters.plotters import plot_form_xz
from compas_thrust.plotters.plotters import plot_form_joints
from numpy import array

from compas.geometry import intersection_segment_segment_xy
from compas.geometry import is_intersection_segment_segment_xy

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file_save = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/hor_action.json'
    lambda_hor = 0.0
    exitflag = 0
    blocks = 6
    form = create_arch(total_nodes = blocks, lambda_hor = lambda_hor, D = 2.00)
    print_opt = True
    thk = 0.25
    form = circular_joints(form, thk = thk, blocks = blocks)
    # plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()

    # Initial parameters

    translation = form.attributes['tmax']
    bounds_width = 5.0
    use_bounds = False
    qmax = 20000
    indset = None

    # Optimisation

    exitflag = 1
    count = 0

    while count < 100 and exitflag is not 0:
        # form = _form(form, keep_q=True)
        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver='slsqp',
                                            printout=print_opt,
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
            plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()
            # form.to_json(file_save)

        i_k = form.index_key()
        joints = form.attributes['joints']
        result = []
        for i, elements in joints.items():
            p1, p2, edges = elements
            # print('Joint: ',i)
            joint_int = False
            for k1, k2 in edges:
                e1 = form.vertex_coordinates(i_k[k1])[0], form.vertex_coordinates(i_k[k1])[2]
                e2 = form.vertex_coordinates(i_k[k2])[0], form.vertex_coordinates(i_k[k2])[2]
                edge_2D = [e1,e2]
                joint_2D = [[p1[0],p1[2]],[p2[0],p2[2]]]
                # print('Edge',k1,k2)
                if is_intersection_segment_segment_xy(joint_2D,edge_2D) == True:
                    joint_int = True
                    break
            result.append(joint_int)

        print(result)
        
        plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True, joints=True).show()
    
    lambda_hor = lambda_hor + 0.02
