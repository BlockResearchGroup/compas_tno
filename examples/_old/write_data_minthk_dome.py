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

    # type_vault = 'cross'
    # type_fd = 'fan_fd'
    # example = 'rectangular/7,5x10/'

    type_vault = 'dome'
    type_fd = 'radial_spaced'
    example = 'xxxx'
    objective = 'max'
    thck = 0.50
    decrease = 0.01
    min_t = 0.207
    R = 5.0
    n_radial = 8
    n_spikes = 20
    divisions = 16

    if type_vault == 'dome':
        FILECSV = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + '/' + type_fd + '/' + 'minthck_via_' + objective + '_' + type_fd + '.csv'
    if type_vault in ['cross','fan']:
        FILECSV = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + '/' + example + type_fd + '/' + 'minthck_via_' + objective + '_' + type_fd + '.csv'

    with open(FILECSV, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Thickness", "fopt/W", "t/R", "Qmax", "Zb", "fopt", "output"])

        while thck > min_t:

            if type_vault == 'dome':
                PATH = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + '/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)
                # file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'
            if type_vault in ['cross','fan']:
                PATH = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + '/' + example + type_fd + '/' + type_fd + '_discr_'+ str(divisions)
                # file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'

            file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'

            print('---------------------- Retrieving Results t = {0}'.format(thck))
            print(file_open)
            form = FormDiagram.from_json(file_open)
            pzt = 0.0
            for key in form.vertices():
                pzt +=form.vertex_attribute(key, 'pz')
            for key in form.vertices_where({'is_fixed': True}):
                rx = round(form.vertex_attribute(key, '_rx'),3)
                ry = round(form.vertex_attribute(key, '_ry'),3)
                zb = round(form.vertex_attribute(key,'z'),3)
                break
            exitflag = form.attributes['exitflag']
            q = [form.edge_attribute(key, 'q') for key in form.edges()]
            fopt = round(form.attributes['fopt'],3)
            qmax = round(max(q),3)
            f_adim = round(form.attributes['fopt']/pzt,4)

            writer.writerow([str(int(round(thck*100))), f_adim, round(thck/R,4), qmax, zb, fopt, exitflag])

            thck = round(thck - decrease, 3)

        thck = min_t
        file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*1000))) + '.json'


        print('---------------------- Retrieving Results t = {0}'.format(thck))
        print(file_open)
        form = FormDiagram.from_json(file_open)
        pzt = 0.0
        for key in form.vertices():
            pzt +=form.vertex_attribute(key, 'pz')
        for key in form.vertices_where({'is_fixed': True}):
            rx = round(form.vertex_attribute(key, '_rx'),3)
            ry = round(form.vertex_attribute(key, '_ry'),3)
            zb = round(form.vertex_attribute(key,'z'),3)
            break
        exitflag = form.attributes['exitflag']
        q = [form.edge_attribute(key, 'q') for key in form.edges()]
        fopt = round(form.attributes['fopt'],3)
        qmax = round(max(q),3)
        f_adim = round(form.attributes['fopt']/pzt,4)

        writer.writerow([str(int(round(thck*100))), f_adim, round(thck/R,4), qmax, zb, fopt, exitflag])

        print('saved on:', FILECSV)
