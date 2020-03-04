from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
# from compas_tno.algorithms.equilibrium import reactions
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_force
from compas_tno.algorithms import initialize_problem
from compas_tno.algorithms import update_tna
from compas_tno.algorithms import z_from_form
from compas_tno.diagrams.form import overview_forces
from compas.utilities import geometric_key

from compas_tno.utilities import fix_boundaries_sym
from compas_tno.utilities import fix_boundaries_complete

from compas_tno.utilities import fix_mid_sym
from compas_tno.utilities import fix_mid_complete

from compas.geometry import is_point_on_segment
from compas.geometry import intersection_segment_segment
from compas.geometry import is_intersection_line_line
from compas.geometry import intersection_line_line
from numpy import shape
from numpy import abs
from numpy import argmin
from numpy import array
from numpy import float64
from numpy import dot
from numpy import hstack
from numpy import isnan
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum
from numpy import vstack
from numpy import zeros
from numpy import ones
from numpy import append
from numpy.linalg import pinv
from numpy.linalg import matrix_rank
from numpy.random import rand
from numpy.random import randint

from compas.numerical import normrow
# from compas.numerical import norm

from scipy.linalg import svd
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

import matplotlib
import matplotlib.pyplot as plt


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    i_s = []
    n_s = []
    lp_s = []

    for i in range(2,12):
        j = 1

        file_load_sym = '/Users/mricardo/compas_dev/me/discretize/0'+ str(j) +'_0'+ str(i) +'_sym.json'
        form = FormDiagram.from_json(file_load_sym)
        form = fix_mid_sym(form, plot=True)
        file_save_sym = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+ str(j) +'_0'+ str(i) +'_sym.json'
        form.to_json(file_save_sym)

        file_load_complete = '/Users/mricardo/compas_dev/me/discretize/0'+ str(j) +'_0'+ str(i) +'_complete.json'
        form = FormDiagram.from_json(file_load_complete)
        form = fix_mid_complete(form, plot=True)
        file_save_complete = '/Users/mricardo/compas_dev/me/loadpath/midsupport/discretize/0'+ str(j) +'_0'+ str(i) +'_complete.json'
        form.to_json(file_save_complete)
