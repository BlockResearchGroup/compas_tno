from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters.plotters import plot_form

import math
import csv


def import_csv(file_csv):

    xmin = []
    xmax = []
    fmin = []
    fmax = []

    with open(file_csv) as data:
        csv_reader = csv.reader(data, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                xmin.append(float(row[0]))
                xmax.append(float(row[2]))
                fmin.append(float(row[1])*100)
                fmax.append(float(row[3])*100)
                line_count += 1
        print(f'Processed {line_count} lines.')


    return xmin, xmax, fmin, fmax


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file_csv1 = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/cross_fd.csv'
    file_csv2 = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/fan_fd.csv'
    # file_csv3
    file_save = '/Users/mricardo/Documents/ETH/Conferences/SAHC2020/My Paper/pictures/Figure_3.pdf'

    xmin1, xmax1, fmin1, fmax1 = import_csv(file_csv1)
    xmin2, xmax2, fmin2, fmax2 = import_csv(file_csv2)
    # xmin3, xmax3, fmin3, fmax3 = import_csv(file_csv3)
    print(xmin1)
    import matplotlib.pyplot as plt
    from numpy import arange
    from matplotlib.ticker import FormatStrFormatter

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    max_x = 0.10
    min_x = 0.05
    interval_x = 0.01
    max_y = 100
    min_y = 60
    interval_y = 10

    x1_ = xmin1[len(xmin1)-1]
    y1_ = fmin1[len(xmin1)-1]
    x2_ = xmin2[len(xmin2)-1]
    y2_ = fmin2[len(xmin2)-1]

    print(x1_, y1_, x2_, y2_)

    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmin1, fmin1, 'o', ls='-', markersize=5, color='blue', label='minimum thrust\n(cross)')
    ax.plot(xmax1, fmax1, 'o', ls='-', markersize=5, color='red', label='maximum thrust\n(cross)')

    ax.plot(xmin2, fmin2, 'x', ls='-', markersize=6, color='blue', label='minimum thrust\n(fan)') # '#009900'
    ax.plot(xmax2, fmax2, 'x', ls='-', markersize=6, color='red', label='maximum thrust\n(fan)') # '#FF8000'

    ax.plot(x1_, y1_, 'o', ls=' ', markersize=7, color='black', label='limit state\n(cross)')
    ax.plot(x2_, y2_, 'X', ls=' ', markersize=8, color='black', label='limit state\n(fan)')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([      0,     33.33,   57.14,    75.00,   88.9,  100])
    ax2.set_xticklabels(['1.0',   '1.2',  '1.4',    '1.6',  '1.8',  '2.0'], size=size_axis_data)
    ax1.set_xlabel('t/R', size=size_axis_label, weight='bold', labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('T/W (%)', size=size_axis_label, weight='bold', labelpad=8)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

    # plt.axvline(x=(x_ - max_x)/(min_x - max_x)*100, ls='--', color='black') # add    ymin=(max_y-y_)/(max_y-min_y)

    box = ax.get_position()
    ax.set_position([box.x0*0.5, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1.03, 0.5), fontsize=size_legend)    #, bbox_to_anchor(0.1, 0.1), ncol=1)

    plt.savefig(file_save, transparent=True)
    plt.show()

