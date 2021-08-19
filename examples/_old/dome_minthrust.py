from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import create_dome_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights
from compas_tno.utilities.constraints import set_dome_heights

from compas.utilities import geometric_key

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters import plot_form
from compas_tno.utilities import check_constraints

from compas_tno.algorithms import z_update

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    for px in [0.35]:

        # Try with 'radial' and 'flower' and for the objective change 'min' and 'max'
        type_fd = 'par-diag'
        objective = 'max'
        thck = 0.30
        px = - 1* px

        # Create Vault from one of the patterns Fan/Grid with the dimensions or load in case of flower FD

        xc = 5.0
        yc = 5.0
        radius = 5.0
        n_radial = 8
        n_spikes = 16

        if type_fd != 'radial':
            PATH = '/Users/mricardo/compas_dev/me/minmax/dome/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes) + '_px_' + str(px)
            # form = FormDiagram.from_json(PATH + '.json')
            # print('Loaded: ', PATH, '.json')
            if type_fd == 'diag':
                form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0, diagonal=True)
            if type_fd == 'par-diag':
                form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0, diagonal=True, partial=True)
        else:
            PATH = '/Users/mricardo/compas_dev/me/minmax/dome/r=' + str(int(radius)) + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)
            form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)

        form = set_dome_heights(form, center=[xc, yc], radius=radius, thck=thck)
        form = delete_boundary_edges(form)

        # plot_form(form).show()

        # Copy Correct Loads:

        # tchk_radial = thck
        # type_fd = 'radial'
        # PATH_rad = '/Users/mricardo/compas_dev/me/minmax/dome/r=' + str(int(radius)) + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)
        # file_radial = PATH_rad + '_' + objective + '_t=' + str(int(round(tchk_radial*100))) + '.json'
        # print(file_radial)

        PATHspecial = '/Users/mricardo/compas_dev/me/minmax/dome/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes) + '_px_' + str(-0.35)
        file_initial_special = PATHspecial + '_' + 'min' + '_t=' + str(int(round(thck*100))) + '.json'

        form_radial = FormDiagram.from_json(file_initial_special)

        gkey_pz = {}
        mid_q = {}
        gkey_zb = {}

        for key in form_radial.vertices():
            pz_radial = form_radial.vertex_attribute(key, 'pz')
            gkey = geometric_key(form_radial.vertex_coordinates(key)[:2] + [0])
            gkey_pz[gkey] = pz_radial
            if form_radial.vertex_attribute(key, 'is_fixed') == True:
                gkey_zb[gkey] = form_radial.vertex_attribute(key, 'z')

        for u,v in form_radial.edges():
            q_radial = form_radial.edge_attribute((u,v), 'q')
            gkey = geometric_key(form_radial.edge_midpoint(u,v)[:2] + [0])
            mid_q[gkey] = q_radial

        for key in form.vertices():
            gkey = geometric_key(form.vertex_coordinates(key)[:2] + [0])
            try:
                pz_radial = gkey_pz[gkey]
                form.vertex_attribute(key, 'pz', pz_radial)
                form.vertex_attribute(key, 'px', pz_radial*px)
            except BaseException:
                print(gkey)
            if form.vertex_attribute(key, 'is_fixed') == True:
                form.vertex_attribute(key, 'z', gkey_zb[gkey])

        for u,v in form.edges():
            gkey = geometric_key(form.edge_midpoint(u,v)[:2] + [0])
            try:
                q_radial = mid_q[gkey]
                form.edge_attribute((u,v), 'q', q_radial)
            except BaseException:
                form.edge_attribute((u,v), 'q', 0.0)

        file_initial = PATH + '_lp.json'

        # objective_load = 'min'
        # file_initial = PATH + '_' + objective_load + '_t=' + str(int(round(thck*100))) + '.json'

        # temp to start from minimum

        # file_initial = PATH + '_' + 'min' + '_t=' + str(int(thck*100)) + '.json'

        file_save = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'

        translation = True
        qmax = 100
        qmin = -1e-6
        indset = None
        print_opt = True

        plot_form(form, show_q=False, fix_width=False).show()

        # plot_form(form, show_q = False, fix_width= True).show()

        # Convex Optimisation to find good starting point. Save the starting point, and can load it later if wanted

        # fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
        #                                             printout=print_opt,
        #                                             find_inds=True, plot = False,
        #                                             tol=0.01,
        #                                             objective='loadpath',
        #                                             indset=indset)
        # form.to_json(file_initial)
        # form = FormDiagram.from_json(file_initial)
        # form = z_update(form)
        # check_constraints(form, show = True)
        # indset = form.attributes['indset']
        # plot_form(form, show_q=False, fix_width=False).show()

        # Maximum or Minimum Thrusts

        solver = 'pyOpt-' + 'SLSQP'
        # solver = 'slsqp'

        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                                    printout=print_opt,
                                                    find_inds=True,
                                                    indset=indset,
                                                    translation=translation,
                                                    objective=objective,
                                                    bmax=False,
                                                    summary=print_opt)

        form.to_json(file_save)
        overview_forces(form)
        # plot_form(form, show_q=False).show()
        reactions(form)

        for key in form.vertices_where({'is_fixed': True}):
            rx = round(form.vertex_attribute(key, '_rx'), 3)
            ry = round(form.vertex_attribute(key, '_ry'), 3)
            zb = round(form.vertex_attribute(key, 'z'), 3)
            print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
            break
        qmax = round(max(qopt).item(), 3)
        fopt = round(fopt, 3)
        print('zb: {0:.3f}'.format(zb))
        print('fopt: {0:.3f}'.format(fopt))
        print('qmax: {0:.3f}'.format(max(qopt).item()))
