from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form

from compas_tno.utilities.constraints import set_pavillion_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters import plot_form

import math

import csv


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'cross_fd'
    objective = 'min'
    total = 10

    thck = 0.30
    reduction = 0.01

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

    # Open initial formdiagram and output file

    PATH = '/Users/mricardo/compas_dev/me/minmax/pavillion/'+ type_fd + '/' + type_fd + '_discr_'+ str(divisions)
    FILECSV = '/Users/mricardo/compas_dev/me/minmax/pavillion/minthck_via_' + objective + '_' + type_fd + '.csv'

    file_initial = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'
    form = FormDiagram.from_json(file_initial)
    for key in form.vertices_where({'is_fixed': True}):
        rx = round(form.vertex_attribute(key, 'rx'),3)
        ry = round(form.vertex_attribute(key, 'ry'),3)
        zb = round(form.vertex_attribute(key,'z'),3)
        break
    q = [form.edge_attribute(key, 'q') for key in form.edges()]
    fopt = form.attributes['fopt']

    with open(FILECSV, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Thickness", "Reaction X", "Reaction Y", "Qmax", "Zb", "fopt", "output"])
        writer.writerow([thck*100, rx, ry, max(q), zb, fopt, 1.0])

    # Convex Optimisation to find good starting point. Save the starting point, and can load it later if wanted

        exitflag = 1
        i = 1
        while exitflag == 1 and i <= total:

            file_initial = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'
            form = FormDiagram.from_json(file_initial)

            thck = round(thck - reduction, 3)
            print('----------------------\nOptimisation with thickness: {0}'.format(thck))
            file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'
            form = set_pavillion_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], thk = thck, b = 5.0, set_heights=False, ub_lb = True, update_loads = True)
            indset = form.attributes['indset']

            # Initial parameters

            translation = form.attributes['tmax']
            bounds_width = 5.0
            use_bounds = False
            qmax = 100
            qmin = -1e-6
            indset = form.attributes['indset']
            print_opt = False

            # Maximum or Minimum Thrusts

            solver = 'pyOpt-' + 'SLSQP'

            fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                                printout=print_opt,
                                                find_inds=True,
                                                indset=indset,
                                                tol=0.01,
                                                translation = translation,
                                                objective=objective,
                                                bmax = True,
                                                summary=print_opt)

            print('fopt: {0:.3f}'.format(fopt))
            overview_forces(form)
            # plot_form(form, show_q = False).show()
            i += 1

            if exitflag == 1:
                form.to_json(file_save)
                for key in form.vertices_where({'is_fixed': True}):
                    rx = round(form.vertex_attribute(key, 'rx'),3)
                    ry = round(form.vertex_attribute(key, 'ry'),3)
                    zb = round(form.vertex_attribute(key,'z'),3)
                    # print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
                    # print('zb: {0:.3f}'.format(zb))
                    break
                qmax = round(max(qopt).item(),3)
                fopt = round(fopt,3)
                # print('fopt: {0:.3f}'.format(fopt))
                # print('qmax: {0:.3f}'.format(max(qopt).item()))
                # exitflag = 0

                writer.writerow([thck*100, rx, ry, qmax, zb, fopt, exitflag])
