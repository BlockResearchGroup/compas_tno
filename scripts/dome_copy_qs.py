from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import create_dome_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights
from compas_tno.utilities.constraints import set_dome_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters import plot_form

from compas_tno.utilities import check_constraints

from compas.utilities import geometric_key

import math

import csv


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'diag'
    objective = 'min'
    thck = 0.08
    reduction = 0.01
    total = 25

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    xc = 5.0
    yc = 5.0
    radius = 5.0
    n_radial = 8
    n_spikes = 16

    no_diag_file = '/Users/mricardo/compas_dev/me/minmax/dome/r=5/radial_discr_8_16_min_t=8.json'
    old_form = FormDiagram.from_json(no_diag_file)

    form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0, diagonal=True)
    form = set_dome_heights(form, center=[xc, yc], radius=radius, thck=thck)
    form = delete_boundary_edges(form)
    PATH = '/Users/mricardo/compas_dev/me/minmax/dome/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)

    file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'

    translation = True
    qmax = 100
    qmin = -1e-6
    indset = None
    print_opt = True

    # plot_form(old_form, show_q=False, fix_width=False).show()
    # plot_form(form, show_q=False, fix_width=False).show()

    check_constraints(old_form, show=True)

    old_qs = {}
    old_zbs = {}

    for u, v in old_form.edges():
        gkey = geometric_key(old_form.edge_midpoint(u, v)[:2] + [0.0])
        old_qs[gkey] = old_form.edge_attribute((u, v), 'q')
    for key in old_form.vertices_where({'is_fixed': True}):
        gkey = geometric_key(old_form.vertex_coordinates(key)[:2] + [0.0])
        old_zbs[gkey] = old_form.vertex_attribute(key, 'z')

    for u, v in form.edges():
        gkey = geometric_key(form.edge_midpoint(u,v))
        try:
            old_q = old_qs[gkey]
            form.edge_attribute((u,v), 'q', old_q)
        except:
            form.edge_attribute((u,v), 'q', 0.0)
    for key in form.vertices_where({'is_fixed': True}):
        gkey = geometric_key(form.vertex_coordinates(key))
        zb = old_zbs[gkey]
        form.vertex_attribute(key, 'z', zb)

    print(len(old_qs))
    print(old_form.number_of_edges())

    # plot_form(form, show_q=False, fix_width=False).show()

    print(file_save)
    form.to_json(file_save)


