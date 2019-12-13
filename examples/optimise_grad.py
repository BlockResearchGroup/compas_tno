from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.grad_based import optimise_tna
from compas_tno.algorithms.scale import evaluate_scale
from compas_tno.algorithms.scale import lagrangian_scale
from compas_tno.algorithms.scale import scale_form
from compas_tno.algorithms.equilibrium import z_from_form

from compas_tno.utilities.constraints import check_constraints
from compas_tno.utilities.symmetry import replicate

from compas_tno.diagrams.form import energy
from compas_tno.diagrams.form import loadpath
from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import adapt_tna
from compas_tno.diagrams.form import evaluate_a


from compas_tno.plotters.plotters import plot_form

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/me/bestfit/pillow/pillow3_init.json'
    file_save = '/Users/mricardo/compas_dev/me/bestfit/pillow/pillowRV_calc.json'

    form = FormDiagram.from_json(file_save)
    # form = z_from_form(form)
    plot_form(form).show()

    # form = adapt_tna(form, zmax = 5.0, method = 'nodal', plot = False, delete_face = False, kmax = 100)
    # form = z_from_form(form)
    # plot_form(form).show()

    r = evaluate_scale(form, loadpath, [0.5,5.5], plot = False)
    print('Scaling of {0}'.format(r))
    form = scale_form(form,r)
    f_init = loadpath(form)
    print('Initial ENERGY after scaling {0}'.format(f_init))
    evaluate_a(form)

    form = optimise_tna(form, objective='loadpath', plot=False, it_max=500, alpha=1.0, a_max = 2.0, steplength=None, null=None, save_steps=False)

    plot_form(form).show()
    f = loadpath(form)
    print('Final energy {0}'.format(f))

    if f < f_init:
        print('Optimisation improved! Var: {0:.3f}%'.format(100*(f-f_init)/f_init))
        form.to_json(file_save)
    else:
        print('It can\'t get better than this!')


