from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces

from compas_tno.algorithms.equilibrium import reactions

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
    thck = 0.50
    reduction = 0.01
    minimum = 0.17

    # Change to 'pavillion/' to 'cross/square' or 'rectangular/5x10' or 'rectangular/7,5x10'
    type_vault = 'pavillion/'

    if type_fd == 'cross_fd':
        divisions = 20
    if type_fd == 'fan_fd':
        divisions = 1

    # Open initial formdiagram and output file

    PATH = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + type_fd + '/' + type_fd + '_discr_'+ str(divisions)
    FILECSV = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + 'minthck_via_' + objective + '_' + type_fd + '.csv'

    print('----------------------\nOptimisation with thickness: {0}'.format(thck))

    file_initial = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'
    form = FormDiagram.from_json(file_initial)

    for key in form.vertices_where({'is_fixed': True}):
        rx = round(form.vertex_attribute(key, 'rx'),3)
        ry = round(form.vertex_attribute(key, 'ry'),3)
        zb = round(form.vertex_attribute(key,'z'),3)
        break
    exitflag = form.attributes['exitflag']
    q = [form.edge_attribute(key, 'q') for key in form.edges()]
    fopt = round(form.attributes['fopt'],3)
    qmax = round(max(q),3)

    with open(FILECSV, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Thickness", "Reaction X", "Reaction Y", "Qmax", "Zb", "fopt", "output"])
        writer.writerow([thck*100, rx, ry, qmax, zb, fopt, exitflag])

        while thck > minimum:

            thck = round(thck - reduction, 3)
            file = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'
            form = FormDiagram.from_json(file)

            print('----------------------\nOptimisation with thickness: {0}'.format(thck))

            for key in form.vertices_where({'is_fixed': True}):
                rx = round(form.vertex_attribute(key, 'rx'),3)
                ry = round(form.vertex_attribute(key, 'ry'),3)
                zb = round(form.vertex_attribute(key,'z'),3)
                break

            q = [form.edge_attribute(key, 'q') for key in form.edges()]
            fopt = form.attributes['fopt']
            exitflag = form.attributes['exitflag']
            qmax = round(max(q),3)
            fopt = round(fopt,3)
            writer.writerow([thck*100, rx, ry, qmax, zb, fopt, exitflag])
