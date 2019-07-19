from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.grad_based import optimise_tna
from compas_thrust.algorithms.grad_based import evaluate_scale
from compas_thrust.algorithms.grad_based import lagrangian_scale
from compas_thrust.algorithms.grad_based import energy

from compas_thrust.algorithms.equilibrium import scale_form

from compas_thrust.utilities.utilities import check_constraints
from compas_thrust.utilities.utilities import oveview_forces
from compas_thrust.utilities.utilities import replicate

from compas_thrust.diagrams.form import adapt_tna

from compas_thrust.plotters.plotters import plot_form

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/me/bestfit/pillowRV_grad.json'
    file_save = '/Users/mricardo/compas_dev/me/bestfit/pillowRV_grad.json'

    form = FormDiagram.from_json(file)

    form = adapt_tna(form, zmax = 5.0, plot = False, delete_face = True)

    r = evaluate_scale(form, energy, [0.5,5.5], plot = True)
    print('Scaling of {0}'.format(r))
    form = scale_form(form,r)
    f_init = energy(form)
    print('Initial ENERGY after scaling {0}'.format(f_init))

    form = optimise_tna(form, plot=False, it_max=5, alpha=1.0, a_max = 2.0, steplength=None, null=None, save_steps=False)

    plot_form(form).show()
    f = energy(form)
    print('Final energy {0}'.format(f))
    
    if f < f_init:
        print('Optimisation improved! Var: {0:.3f}%'.format(100*(f-f_init)/f_init))
        form.to_json(file_save)
    else:
        print('It can\'t get better than this!')


