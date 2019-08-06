
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single
from compas_thrust.algorithms.ind_based import initialize_problem

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.diagrams.form import oveview_forces
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import _form
from compas_thrust.plotters.plotters import plot_form

from compas_thrust.plotters.plotters import plot_form
import compas_pattern

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":


    file = '/Users/mricardo/compas_dev/me/minmax/radial/02_07_complete.json'
    # file_save = '/Users/mricardo/compas_dev/me/minmax/2D_arch/01_joints_min.json'
    # file_complete = '/Users/mricardo/compas_dev/me/minmax/radial/02_06_complete.json'

    form = FormDiagram.from_json(file)

    # Initial parameters

    # tmax = 0.5 # form.attributes['tmax']
    # bounds_width = 2.5
    # use_bounds = False
    # qmax = 120
    # indset = None

    plot_form(form,radius=0.005).show()

    # # Optimisation

    # fopt, qopt = optimise_single(form, qmax=qmax, solver='devo',
    #                                     polish='slsqp',
    #                                     population=100,
    #                                     generations=50,
    #                                     printout=50,
    #                                     tol=0.01,
    #                                     t = tmax,
    #                                     opt_max=False,
    #                                     tension=False,
    #                                     use_bounds = use_bounds,
    #                                     bounds_width = bounds_width,
    #                                     objective='min',
    #                                     indset=indset,
    #                                     buttress = True)

    # # Check compression and Save

    # q = [attr['q'] for u, v, attr in form.edges(True)]
    # qmin  = min(array(q))
    # if qmin > -0.1: # and check_constraints(form) < 1.0:
    #     oveview_forces(form)
    #     reactions(form, plot=False)
    #     print('Optimisation completed')
    #     form.to_json(file_save)

    # Replicate-sym and Print Results

    # print('Horizontal checks: {0}'.format(horizontal_check(form)))

    # form_ = replicate(form, file_complete)
    reactions(form)
    check_constraints(form)
    # form.to_json(file)
    # oveview_forces(form_)
    # plot_form(form_)
