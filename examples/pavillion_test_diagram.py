from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters import plot_form

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'cross_fd'
    objective = 'max'
    thck = 0.17

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    x_span = 10.0
    y_span = 10.0

    if type_fd == 'cross_fd':
        divisions = 20
        xy_span = [[0.0,x_span],[0.0,y_span]]
        form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions, fix='all') # FIX ALL NODES ON BOUNDARIES
    if type_fd == 'fan_fd':
        divisions = 16
        xy_span = [[0.0,x_span],[0.0,y_span]]
        form = create_fan_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions, fix='all') # FIX ALL NODES ON BOUNDARIES
    plot_form(form, show_q=False).show()
    form = delete_boundary_edges(form)
    form = set_pavillion_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], thk = thck, b = 5.0, t = 0.0, set_heights=False, ub_lb = True, update_loads = True)
    plot_form(form, show_q=False).show()




