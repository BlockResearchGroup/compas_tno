
from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.ind_based import optimise_single
from compas_tno.algorithms.ind_based import initialize_problem

from compas_tno.algorithms.equilibrium import reactions
from compas_tno.algorithms.equilibrium import horizontal_check

from compas_tno.algorithms import optimise_general

from compas_tno.utilities.constraints import circular_heights
from compas_tno.utilities.constraints import check_constraints
from compas_tno.diagrams.form import overview_forces
from compas_tno.utilities.symmetry import replicate

from compas_viewers.meshviewer import MeshViewer

from compas_tno.diagrams.form import _form
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_form_xz

from compas_tno.plotters import plot_form
import compas_pattern

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":


    # file = '/Users/mricardo/compas_dev/me/minmax/radial/01_05_sym.json'
    file = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/01_lp.json'
    # file_save = '/Users/mricardo/compas_dev/me/minmax/fan/fill_01_05_calc.json'
    # file_complete = '/Users/mricardo/compas_dev/me/minmax/fan/fill_01_05_complete.json'

    # file = '/Users/mricardo/compas_dev/me/convex/4bars/diagram.json'
    px = 1.0
    form = FormDiagram.from_json(file)
    # check_constraints(form, show=True, lb_show=True, ub_show=True)
    thk = 0.25
    exitflag = 0
    while exitflag == 0: # in [0.20]:
        form = circular_heights(form, thk = thk)
        # form = _form(form)
        lb = []
        ub = []
        xs = []
        px = px+ 0.1
        print('\n\n\n-----------------------PX ',px,'\n\n\n')
        for key in form.vertices():
            if key == 9:
                lb.append(form.vertex_attribute(key,'lb'))
                ub.append(form.vertex_attribute(key,'ub'))
                x, y, z = form.vertex_coordinates(key)
                form.vertex_attribute(key, 'px', value = px)
            # if form.vertex_attribute(key, 'b') is not None:
                # form.vertex_attribute(key, 'b', value = [0.10, 0.0] )
                # print(form.vertex_attribute(key, 'b'))
            # print(x)
            # except BaseException:
            #     pass
            # xs.append(x)
        #     print('LB: {0} - UB: {1}'.format(form.vertex_attribute(key,'lb'),form.vertex_attribute(key,'ub')))
        # print(len(ub))
        # print(len(lb))
        # print(max(xs))
        # print(min(xs))

        # Initial parameters

        translation = form.attributes['tmax'] # 0.5
        bounds_width = 5.0
        use_bounds = False
        qmax = 200
        indset = None
        # plot_form(form,radius=0.05, heights=True, show_q=False).show()
        # plot_form_xz(form,radius=0.05, heights=True, show_q=False, thk = thk).show()

        # # Optimisation

        exitflag = 1
        count = 0

        while count < 100 and exitflag is not 0:
            # form = _form(form, keep_q=True)
            fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver='slsqp',
                                                printout=False,
                                                find_inds=True,
                                                tol=0.01,
                                                translation = translation,
                                                tension=False,
                                                use_bounds = use_bounds,
                                                bounds_width = bounds_width,
                                                objective='constr_lp',
                                                indset=indset,
                                                bmax = True,
                                                summary=True)

            # Check compression and Save

            q = [attr['q'] for u, v, attr in form.edges(True)]
            qmin  = min(array(q))
            count += 1
            if qmin > -0.1 and exitflag == 0:
                # overview_forces(form)
                # reactions(form, plot=False)
                print('Optimisation completed - Trial:',count, 'PX: ', px)
                plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True).show()
                # form.to_json(file_save)

        # Replicate-sym and Print Results

        # print('Horizontal checks: {0}'.format(horizontal_check(form)))
        # overview_forces(form)

        # form_ = replicate(form, file_complete)
        # reactions(form_)
        # check_constraints(form_, show=True)
        # form.to_json(file)
        # oveview_forces(form_)
        # plot_form_xz(form, radius=0.01, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thk, plot_reactions=True).show()
        # thk = thk - 0.0001
    # viewer = MeshViewer()
    # viewer.mesh = form
    # viewer.show()
