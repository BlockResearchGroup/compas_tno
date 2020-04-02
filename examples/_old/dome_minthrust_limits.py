from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import create_dome_form
from compas_tno.diagrams.form import create_dome_form_spaced
from compas_tno.diagrams.form import create_dome_flower
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

    exitflag = 0
    thck = 0.50
    reduction = 0.01
    data_thk = []
    data_fobj = []
    data_exitflag = []
    i = 1

    while exitflag == 0:
    # while i < 2:

        data_thk.append(thck)

        for objective in ['max']:

            # Try with 'radial' and 'flower' and for the objective change 'min' and 'max'
            type_fd = 'radial'
            # objective = 'max'

            # Create Vault from one of the patterns Fan/Grid with the dimensions or load in case of flower FD

            xc = 5.0
            yc = 5.0
            radius = 5.0
            n_radial = 8
            n_spikes = 20

            PATH = '/Users/mricardo/compas_dev/me/minmax/dome/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)

            if i >= 1:
                thck_load = thck
                objective_load = objective
                if reduction > 0.001:
                    form_load = PATH + '_' + objective_load + '_t=' + str(int(round(thck_load*100))) + '.json'
                else:
                    form_load = PATH + '_' + objective_load + '_t=' + str(int(round(thck_load*1000))) + '.json'
                print('Load:', form_load)
                form = FormDiagram.from_json(form_load)
                check_constraints(form, show=True)
                indset = form.attributes['indset']
                thck = round(thck - reduction, 4)
                print('---------------- THCK: ', thck)
            else:
                print('---------------- THCK: ', thck)
                if type_fd == 'radial':
                    form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)
                elif type_fd == 'flower':
                    form = create_dome_flower(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)
                elif type_fd == 'diag':
                    form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0, diagonal=True)
                elif type_fd == 'par-diag':
                    form = create_dome_form(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0, diagonal=True, partial=True)
                elif type_fd == 'radial_spaced':
                    form = create_dome_form_spaced(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)
            # form = FormDiagram.from_json(PATH + '_lp.json')
            # form = create_dome_flower(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)
            form = set_dome_heights(form, center=[xc, yc], radius=radius, thck=thck)
            form = delete_boundary_edges(form)

            # plot_form(form).show()

            # Copy Correct Loads:

            # file_initial = PATH + '_lp.json'

            # objective_load = 'min'
            # file_initial = PATH + '_' + objective_load + '_t=' + str(int(round(thck*100))) + '.json'

            # temp to start from minimum

            file_initial = PATH + '_lp.json'
            if reduction > 0.001:
                file_save = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'
            else:
                file_save = PATH + '_' + objective + '_t=' + str(int(round(thck*1000))) + '.json'
            print(file_save)

            translation = True
            qmax = 10000
            qmin = -1e-6
            print_opt = True
            bmax = True

            # plot_form(form, show_q=False, fix_width=False).show()

            # Convex Optimisation to find good starting point. Save the starting point, and can load it later if wanted

            if i < 1:
                fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
                                                            printout=print_opt,
                                                            find_inds=True, plot = False,
                                                            tol=0.01,
                                                            objective='loadpath',
                                                            indset=None)
                form.to_json(file_initial)
                indset = form.attributes['indset']
                # plot_form(form, show_q=False, fix_width=False).show()


            for key in form.vertices_where({'is_fixed': True}):
                ub = form.vertex_attribute(key, 'ub')
                form.vertex_attribute(key, 'z', value = ub)

            # Maximum or Minimum Thrusts

            # solver = 'pyOpt-' + 'SLSQP'
            solver = 'slsqp'

            fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                                        printout=print_opt,
                                                        find_inds=True,
                                                        indset=indset,
                                                        translation=translation,
                                                        objective=objective,
                                                        bmax=bmax,
                                                        summary=print_opt)

            reactions(form)
            if exitflag == 0:
                form.to_json(file_save)
            overview_forces(form)
            # plot_form(form, show_q=False).show()
            # reactions(form)

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
            data_fobj.append(fopt)
            data_exitflag.append(exitflag)
            plot_form(form, show_q=False, fix_width=False).show()
        i = i + 1

    print(data_exitflag)
    print(data_fobj)
    print(data_thk)
