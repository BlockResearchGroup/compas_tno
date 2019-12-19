from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
# from compas_tno.algorithms.equilibrium import reactions
from compas_tno.plotters.plotters import plot_form
from compas_tno.plotters.plotters import plot_force
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

    for i in range(5,11):

        print('\n\n------------------ Form ',str(i),'\n')
        j=1

        file_complete = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0'+ str(j) +'_0'+ str(i) +'_complete_paper.json'
        form = FormDiagram.from_json(file_complete)
        # overview_forces(form)
        # i_s.append(i)
        # n_s.append(form.number_of_edges())
        # lp_s.append(form.attributes['loadpath'])
        form = initialize_problem(form)
        plot_form(form,show_q=False, simple=False, max_width=5.0).show()
        force = ForceDiagram.from_formdiagram(form)
        form, force = update_tna(form, delete_face=False)
        plot_force(force, form, color_inds=True).show()

    # fig, ax = plt.subplots()
    # ax.plot(n_s, lp_s)
    # ax.set(xlabel='# Edges F.D.', ylabel='LP',
    #     title='LP for different discretization - Orthogonal')
    # ax.grid()
    # plt.show()

    # print(lp_s)
