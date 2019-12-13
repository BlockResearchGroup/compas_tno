
from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.ind_based import optimise_single
from compas_tno.algorithms.ind_based import initialize_problem

from compas_tno.algorithms.equilibrium import reactions
from compas_tno.algorithms.equilibrium import horizontal_check

from compas_tno.utilities.constraints import check_constraints
from compas_tno.diagrams.form import overview_forces
from compas_tno.utilities.symmetry import replicate

from compas_viewers.meshviewer import MeshViewer

from compas_tno.diagrams.form import _form
from compas_tno.plotters.plotters import plot_form

from compas_tno.plotters.plotters import plot_form
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

    form = FormDiagram.from_json(file)
    # check_constraints(form, show=True, lb_show=True, ub_show=True)

    # form = _form(form)

    # for key in form.vertices():
    #     print('LB: {0} - UB: {1}'.format(form.get_vertex_attribute(key,'lb'),form.get_vertex_attribute(key,'ub')))

    # Initial parameters

    tmax = None
    bounds_width = 5.0
    use_bounds = False
    qmax = 100
    indset = None

    # plot_form(form,radius=0.1, heights=True, show_q=False).show()

    # # Optimisation

    fopt, qopt = optimise_single(form,  qmax=qmax, solver=None,
                                        polish='slsqp',
                                        population=800,
                                        generations=500,
                                        printout=50,
                                        tol=0.01,
                                        t = tmax,
                                        opt_max=False,
                                        tension=False,
                                        use_bounds = use_bounds,
                                        bounds_width = bounds_width,
                                        objective='loadpath',
                                        indset=indset,
                                        buttress = False)

    # Check compression and Save

    q = [attr['q'] for u, v, attr in form.edges(True)]
    qmin  = min(array(q))
    if qmin > -0.1 and check_constraints(form) < 1.0:
        overview_forces(form)
        reactions(form, plot=False)
        print('Optimisation completed')
        # form.to_json(file_save)

    # Replicate-sym and Print Results

    # print('Horizontal checks: {0}'.format(horizontal_check(form)))
    overview_forces(form)

    # form_ = replicate(form, file_complete)
    # reactions(form_)
    # check_constraints(form_, show=True)
    # form.to_json(file)
    # oveview_forces(form_)
    plot_form(form)

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()
