from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.grad_based import optimise_tna
from compas_thrust.algorithms.scale import evaluate_scale
from compas_thrust.algorithms.scale import lagrangian_scale
from compas_thrust.algorithms.scale import scale_form
from compas_thrust.algorithms.equilibrium import z_from_form

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import energy
from compas_thrust.diagrams.form import oveview_forces
from compas_thrust.diagrams.form import adapt_tna
from compas_thrust.diagrams.form import evaluate_a


from compas_thrust.plotters.plotters import plot_form

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/me/minmax/fan/01_05_complete.json'
    # file_save = '/Users/mricardo/compas_dev/me/minmax/fan/01_05_target.json'

    form = FormDiagram.from_json(file)
    form = z_from_form(form)
    plot_form(form).show()

    form = adapt_tna(form, zmax = 5.0, method = 'nodal', plot = True, delete_face = True, kmax = 1000)
    form = z_from_form(form)
    plot_form(form).show()

    r = evaluate_scale(form, energy, [0.5,5.5], plot = True)
    print('Scaling of {0}'.format(r))
    form = scale_form(form,r)
    f_init = energy(form)
    print('Initial ENERGY after scaling {0}'.format(f_init))
    evaluate_a(form)

    form = optimise_tna(form, plot=False, it_max=5, alpha=1.0, a_max = 2.0, steplength=None, null=None, save_steps=False)

    plot_form(form).show()
    f = energy(form)
    print('Final energy {0}'.format(f))
    
    if f < f_init:
        print('Optimisation improved! Var: {0:.3f}%'.format(100*(f-f_init)/f_init))
        form.to_json(file_save)
    else:
        print('It can\'t get better than this!')


