from compas_tna.diagrams import FormDiagram

from compas_tno.algorithms.grad_based import optimise_tna
from compas_tno.algorithms.grad_based import evaluate_scale
from compas_tno.algorithms.grad_based import lagrangian_scale
from compas_tno.algorithms.grad_based import energy
from compas_tno.algorithms.grad_based import loadpath

from compas_tno.algorithms.equilibrium import scale_form

from compas_tno.utilities.utilities import check_constraints
from compas_tno.utilities.utilities import oveview_forces
from compas_tno.utilities.utilities import replicate

from compas_tno.algorithms.equilibrium import z_from_form

from compas_tno.diagrams.form import adapt_tna
from compas_tno.diagrams.form import remove_feet
from compas_tno.diagrams.form import evaluate_a

from compas_tno.plotters import plot_form
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
    file_init = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/C_init.json'

    # form = FormDiagram.from_json(file)
    # form = adapt_tna(form, zmax = 5.0, method = 'geometric', plot = False, delete_face = True, alpha = 100.0, kmax = 1000)
    # print('Solution test')
    # evaluate_a(form)
    # plot_form(form).show()
    # form.to_json(file_tna)

    # Inverse

    form = FormDiagram.from_json(file_tna)
    # form.plot()

    # pzt = 0
    # for key, attr in form.vertices(True):
    #     pzt += attr['pz']
    # print('Load applied before calculation: {0}'.format(pzt))

    form = remove_feet(form, plot = True)
    form.to_json(file_init)




