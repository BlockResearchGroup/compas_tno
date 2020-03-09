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

    file_csv = '/Users/mricardo/compas_dev/me/minmax/2D_arch/diagram_thrust/diagram_thrust_arch.csv'
    file_save = '/Users/mricardo/Documents/ETH/Conferences/SAHC2020/My Paper/pictures/Figure_1.pdf'

    xmin = []
    xmax = []
    fmin = []
    fmax = []
    GSF = []

    import csv

    with open(file_csv) as data:
        csv_reader = csv.reader(data, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                xmin.append(float(row[0]))
                xmax.append(float(row[0]))
                GSF.append(float(row[1]))
                fmin.append(float(row[2])*100)
                fmax.append(float(row[3])*100)
                # print(f'Column names are {", ".join(row)}')
                line_count += 1
        print(f'Processed {line_count} lines.')

    # print(x)
    # print(fmin)
    # print(fmax)
    # print(GSF)

    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from numpy import arange
    from matplotlib.ticker import FormatStrFormatter

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    max_x = 0.2
    min_x = 0.1
    interval_x = 0.01
    max_y = 60
    min_y = 20
    interval_y = 10

    x_ = xmin[len(xmin)-1]
    y_ = fmin[len(xmin)-1]

    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='blue', label='minimum thrust')
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='red', label='maximum thrust')
    ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([0,     18.18,  33.33,  46.15,  57.14,  66.66,  75.0,   82.35,  88.88,  94.74,  100])
    ax2.set_xticklabels(['1.0',   '1.1',  '1.2',  '1.3',  '1.4',  '1.5',  '1.6',  '1.7',  '1.8',  '1.9',  '2.0'], size=size_axis_data)
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
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)    #, bbox_to_anchor(0.1, 0.1), ncol=1)

    plt.savefig(file_save, transparent=True)
    plt.show()

