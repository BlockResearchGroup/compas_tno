import math
import csv

from compas_tna.diagrams import FormDiagram
from compas_tno.algorithms import optimise_general

from compas_tno.utilities.constraints import circular_heights

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_arch

from compas_tno.plotters import plot_form_xz

from numpy import array

import subprocess

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    divisions = 20
    blocks = divisions
    lambda_hor = 0.28
    form = create_arch(total_nodes=blocks, total_self_weight=20, lambda_hor=lambda_hor)
    overview_forces(form)
    print_opt = True
    plot_sols = True
    objective = 'min'
    increase = 0.01
    thk = 0.20
    exitflag = 0
    total = 3
    count = 0

    PATH = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/hor-loads/' + 'Blocks_' + str(blocks)
    FILECSV = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/hor-loads/' + '2D_arch_'+ str(blocks) + '_' + objective + '.csv'

    with open(FILECSV, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Lambda Hor", "Reaction X", "Reaction Y", "Qmax", "Zb", "fopt", "output"])

        images = []

        while exitflag == 0 and count < total:

            lambda_hor = round(lambda_hor + increase, 3)
            form = create_arch(total_nodes=blocks, total_self_weight=20, lambda_hor=lambda_hor)
            form = circular_heights(form, thk=thk)
            file_save = PATH + '_px=' + str(lambda_hor) + '_' + objective + '_t=' + str(int(thk*1000)) + '.json'
            file_img = PATH + '_px=' + str(lambda_hor) + '_' + objective + '_t=' + str(int(thk*1000)) + '.tiff'
            file_img = False
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
                print('Optimisation completed - Trial:', count, 'Lambda', lambda_hor)
                if plot_sols:
                    plot_form_xz(form, radius=0.008, save=file_img, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True,).show()
                    images.append(file_img)

                form.to_json(file_save)
                print('File saved: ', file_save)

                for key in form.vertices_where({'is_fixed': True}):
                    rx = round(form.vertex_attribute(key, 'rx'), 3)
                    ry = round(form.vertex_attribute(key, 'ry'), 3)
                    zb = round(form.vertex_attribute(key, 'z'), 3)
                    print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
                    print('zb: {0:.3f}'.format(zb))
                    break
                qmax = round(max(qopt).item(), 3)
                fopt = round(fopt, 3)
                print('fopt: {0:.3f}'.format(fopt))
                print('qmax: {0:.3f}'.format(max(qopt).item()))

                writer.writerow([lambda_hor*100, rx, ry, qmax, zb, fopt, exitflag])

    print('Optimisation with objective ', objective, ' ended with hor_load: ', lambda_hor)
    lambda_hor = round(lambda_hor - increase,3)
    # print('Plotting last feasible solution with lambda:', lambda_hor)
    # file_print = PATH + '_px=' + str(lambda_hor) + '_' + objective + '_t=' + str(int(thk*1000)) + '.json'
    # form = FormDiagram.from_json(file_print)
    # plot_form_xz(form, radius=0.008, simple=True, fix_width=True, max_width=1.5, heights=True, show_q=False, thk=thk, plot_reactions=True,).show()

    # Make a GIF
    # filepath = PATH + '_' + objective + '_GIF.gif'
    # delay = 50
    # loop = 0
    # command = ['convert', '-delay', '{}'.format(delay), '-loop', '{}'.format(loop), '-layers', 'optimize']
    # subprocess.call(command + images + [filepath])

    # Answer: Why pyOpt does not solve and how to deal with the horizontal loads.
