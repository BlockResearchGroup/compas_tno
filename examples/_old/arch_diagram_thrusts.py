import math
import csv

from compas_tna.diagrams import FormDiagram
from compas_tno.algorithms import optimise_general

from compas_tno.utilities.constraints import circular_heights

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_arch

from compas_tno.plotters import plot_form_xz
from compas_tno.plotters import plot_form

from numpy import array

import subprocess

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    divisions = 20
    blocks = divisions
    form = create_arch(total_nodes=blocks, unit_weight=1.0, lambda_hor=None)
    overview_forces(form)
    print_opt = True
    plot_sols = True
    objective = 'max'
    decrease = 0.01
    thk = 0.1079
    exitflag = 0
    total = 125
    count = 0

    form_weight = circular_heights(form, thk=0.20)
    W = form_weight.attributes['selfweight']
    print('Selfweight', W)
    # plot_form(form_weight).show()

    PATH = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/diagram_thrust/' + 'Blocks_' + str(blocks)
    FILECSV = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/diagram_thrust/' + '2D_arch_'+ str(blocks) + '_fine_' + objective + '.csv'

    # file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '_px=' + str(px) + '.json'

    # for key in form.vertices_where({'is_fixed': True}):
    #     rx = round(form.vertex_attribute(key, '_rx'), 3)
    #     ry = round(form.vertex_attribute(key, '_ry'), 3)
    #     zb = round(form.vertex_attribute(key, 'z'), 3)
    #     break
    # q = [form.edge_attribute(key, 'q') for key in form.edges()]
    # fopt = round(form.attributes['fopt'], 3)
    # exitflag = form.attributes['exitflag']

    with open(FILECSV, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Thickness", "Reaction X", "Reaction Y", "Qmax", "Zb", "fopt", "output"])
        # writer.writerow([thk*100, rx, ry, round(max(q), 3), zb, fopt, exitflag])

        images = []

        while exitflag == 0 and count < total:

            form = circular_heights(form, thk=thk, overwrite_weight=W)
            file_save = PATH + '_' + objective + '_t=' + str(int(thk*1000)) + '.json'
            file_img = PATH + '_' + objective + '_t=' + str(int(thk*1000)) + '.tiff'
            print('Radius Ext/Int: ', form.attributes['Re'], form.attributes['Ri'])

            # Initial parameters

            translation = True
            qmax = 20000
            indset = None

            # Optimisation

            solver = 'slsqp'

            fopt, qopt, zbopt, exitflag = optimise_general(form,
                                                           qmax=qmax,
                                                           solver=solver,
                                                           printout=print_opt,
                                                           find_inds=True,
                                                           translation=translation,
                                                           objective=objective,
                                                           indset=indset,
                                                           bmax=True,
                                                           summary=print_opt
                                                           )

            # Check compression and Save

            q = [attr['q'] for u, v, attr in form.edges(True)]
            qmin = min(array(q))
            count += 1
            if qmin > -0.1 and exitflag == 0:

                print('Optimisation completed - Trial:', count, 't', thk)
                if plot_sols:
                    plot_form_xz(form, radius=0.02, save=file_img, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True,).show()
                    images.append(file_img)
                # form.to_json(file_save)
                print('File saved: ', file_save)

                for key in form.vertices_where({'is_fixed': True}):
                    rx = round(form.vertex_attribute(key, '_rx'), 3)
                    ry = round(form.vertex_attribute(key, '_ry'), 3)
                    zb = round(form.vertex_attribute(key, 'z'), 3)
                    print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
                    print('zb: {0:.3f}'.format(zb))
                    break
                qmax = round(max(qopt).item(), 3)
                fopt_ = round(fopt, 3)
                fopt_adim = round(fopt/W,3)
                print('fopt: {0:.3f}'.format(fopt_))
                print('qmax: {0:.3f}'.format(max(qopt).item()))
                print('FOBJ / W: {0:.3f}'.format(fopt_adim))
                writer.writerow([thk*100, rx, ry, qmax, zb, fopt_adim, exitflag])

                thk = round(thk - decrease, 4)

    print('Optimisation with objective ', objective, ' ended with thickness: ', thk)
    thk = round(thk + decrease, 4)
    print('Plotting last feasible solution with thk:', thk)
    file_print = PATH + '_' + objective + '_t=' + str(int(thk*1000)) + '.json'
    form = FormDiagram.from_json(file_print)
    plot_form_xz(form, radius=0.015, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True,).show()

    # Make a GIF
    # filepath = PATH + '_' + objective + '_GIF.gif'
    # delay = 50
    # loop = 0
    # command = ['convert', '-delay', '{}'.format(delay), '-loop', '{}'.format(loop), '-layers', 'optimize']
    # subprocess.call(command + images + [filepath])

    # Answer: Why pyOpt does not solve and how to deal with the horizontal loads.
