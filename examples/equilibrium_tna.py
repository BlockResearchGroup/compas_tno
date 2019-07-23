from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.grad_based import optimise_tna
from compas_thrust.algorithms.grad_based import evaluate_scale
from compas_thrust.algorithms.grad_based import lagrangian_scale
from compas_thrust.algorithms.grad_based import energy
from compas_thrust.algorithms.grad_based import loadpath

from compas_thrust.algorithms.equilibrium import scale_form

from compas_thrust.utilities.utilities import check_constraints
from compas_thrust.utilities.utilities import oveview_forces
from compas_thrust.utilities.utilities import replicate

from compas_thrust.algorithms.equilibrium import z_from_form

from compas_thrust.diagrams.form import adapt_tna
from compas_thrust.diagrams.form import remove_feet
from compas_thrust.diagrams.form import evaluate_a

from compas_thrust.plotters.plotters import plot_form
from compas.utilities import geometric_key

from copy import deepcopy
from numpy import array
from numpy import argmin

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/C_comp.json'
    file_tna = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/C_tna.json'

    form = FormDiagram.from_json(file)
    form = adapt_tna(form, zmax = 5.0, plot = False, delete_face = True, alpha = 0.9)
    print('Solution test')
    evaluate_a(form)

    f_init = loadpath(form)
    print('Initial Loadpath after scaling {0}'.format(f_init))

    # form = optimise_tna(form, plot=False, it_max=5, alpha=1.0, a_max = 2.0, steplength=None, null=None, save_steps=False)

    # plot_form(form).show()
    # form.to_json('/Users/mricardo/compas_dev/compas_loadpath/data/freeform/A_tna.json')

    # Inverse

    # form.plot()
    form = FormDiagram.to_json(file_tna)
    form = remove_feet(form, plot = True)

    


