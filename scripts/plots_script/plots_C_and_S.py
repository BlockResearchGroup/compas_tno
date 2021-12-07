import os
import matplotlib.pyplot as plt
import numpy as np
from compas_tno.plotters import open_csv_row
# from compas_tno.plotters import diagram_of_thrust
# from compas_tno.plotters import diagram_of_multiple_thrust
from matplotlib.ticker import FormatStrFormatter
from numpy import arange
from numpy import array
from numpy import append
from numpy import linspace
from numpy import concatenate
import csv

from compas.geometry import intersection_line_line_xy

# If want to change font
# plt.rcParams["font.family"] = "Times"


def diagram_of_thrust(thicknesses, solutions, limit_state=True, fill=False, xy_limits=None, x_label=None, GSF_ticks=False, save=False):
    """ Plot a diagram of Thrusts based on the collected data from (n) points.

    Parameters
    ----------
    thicknesses : list of lists [[n],[m]]
        Points with discretised solutions for minimum/maximum thrust.
    solutions : list of lists [[n],[m]]
        Adimensional thrust over weight for minimum/maximum thrust.
    limit_state : bool
        If yes, it interpolates the limit_state based on the two last solutions.
    fill : bool
        If yes, it fills bewteen the curves of min and max.
    xy_limits : list of lists
        Freezes the axis of the graphs.
    save : str
        Path to save the graph.

    Returns
    -------
    obj
        Plotter object.

    """
    thicknesses_min, thicknesses_max = thicknesses
    min_sol, max_sol = solutions
    xmin = thicknesses_min
    xmax = thicknesses_max
    fmin = 100.0 * array(min_sol)
    fmax = -100.0 * array(max_sol)
    # print('\n', xmin, xmax, fmin, fmax)
    n = len(xmin)
    m = len(xmax)
    if m > n:
        dimension = xmax
    else:
        dimension = xmin

    if xmin[-1] == xmax[-1]:
        x_ = xmin[-1]
        y_ = fmin[-1]
        limit_state = True
        dim = dimension
        xmax_ = xmax
        xmin_ = xmin
        fmin_ = fmin
        fmax_ = fmax
    else:
        if n >= 2 and m >= 2:
            x_, y_, _ = intersection_line_line_xy([[xmax[m-1], fmax[m-1]], [xmax[m-2], fmax[m-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
        dim = append(dimension, x_)
        xmax_ = append(xmax, x_)
        xmin_ = append(xmin, x_)
        fmin_ = append(fmin, y_)
        fmax_ = append(fmax, y_)

    # print('\n', xmin_, xmax_, fmin_, fmax_)

    # size_axis_label = 14
    # size_axis_data = 12
    # size_legend = 14
    size_axis_label = 12
    size_axis_data = 12
    size_legend = 12

    interval_x = abs(dimension[0] - dimension[1])
    interval_y = 10
    if xy_limits:
        [[max_x, min_x], [max_y, min_y]] = xy_limits
    else:
        max_x = max(dim)
        min_x = min(dim) - interval_x
        max_y = interval_y - max(fmax_) % 10 + max(fmax_)
        min_y = min(fmin_) - min(fmin_) % 10

    if not GSF_ticks:
        last_scale = round(max_x/min_x, 2)
        middle_scale = round((last_scale - 1.0)/2 + 1, 2)
        quad_scale = round((last_scale - 1.0)*3/4 + 1, 2)
        ticks_GSF = [middle_scale, quad_scale, last_scale]
    else:
        ticks_GSF = GSF_ticks

    ticks_x = []
    for tck in ticks_GSF:
        ticks_x.append(max_x/tck)

    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='red', label='maximum thrust')
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='blue', label='minimum thrust')

    if limit_state:
        extrapolation_max = [fmax[m-1], y_]
        extrapolation_min = [fmin[n-1], y_]
        extrapolation_x_max = [xmax[m-1], x_]
        extrapolation_x_min = [xmin[n-1], x_]
        ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')
        ax.plot(extrapolation_x_max, extrapolation_max, '', ls='--', color='red')
        ax.plot(extrapolation_x_min, extrapolation_min, '', ls='--', color='blue')
        ax.annotate(str(round(max_x/x_, 2)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    if fill:
        ax.fill(append(xmin_, xmax_[::-1]), append(fmin_, fmax_[::-1]), color="grey", alpha=0.2)

    ax2 = ax.twiny()

    # ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    # ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF], size=size_axis_data)
    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)

    if x_label:
        ax.set_xlabel(x_label, size=size_axis_label, labelpad=8)
    else:
        ax.set_xlabel(r'thickness ($t$)', size=size_axis_label, labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)
    ax.set_ylabel(r'thrust-over-weight ($T_i/W_i$) [%]', size=size_axis_label, labelpad=8)
    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # Check if this is necessary
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax.tick_params(axis='both', which='major', labelsize=size_axis_data)

    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)
    ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt

def diagram_of_multiple_thrust(thicknesses, solutions, legends, simplified=True, limit_state=True, colors=None, xy_limits=None, GSF_ticks=False, fill=False, show_legend=True, save=None, markers=None):
    """ Plot a diagram of Thrusts based on the collected data on (m) problems each with (n) points.

    Parameters
    ----------
    dimensions : list of lists (m-n)
        Points with discretised solutions.
    min_sol : list of lists (m-n)
        Adimensional thrust over weight for minimum thrust.
    max_sol : list of lists (m-n)
        Adimensional thrust over weight for minimum thrust.
    legends : list (m)
        Legend of the (m) problems.

    Returns
    -------
    obj
        Plotter object.

    """

    kmax = len(legends)

    if markers is None:
        markers = ['o', '^', 's', 'D', 'x', '1', '2', '3', 'v', 'p', '4', '8']

    if colors is None:
        colormap = plt.cm.coolwarm  # gist_ncar nipy_spectral, Set1, Paired coolwarm
        colors = [colormap(i) for i in linspace(0, 1.0, kmax)]

    size_axis_label = 12
    size_axis_data = 12
    size_legend = 12

    # fig = plt.figure(figsize=[15, 5])
    fig = plt.figure(figsize=[12, 5])
    ax = fig.add_subplot(1, 1, 1)

    interval_x = 0.05 # abs(flatten_list(min_thks)[0] - flatten_list(min_thks)[1])
    interval_y = 20

    if xy_limits is None:
        def flatten_list(l): return [item for sublist in l for item in sublist]
        extreme_max = -100*min(flatten_list(max_sols))
        extreme_min = 100*min(flatten_list(min_sols))
        max_x = max(flatten_list(min_thks))
        min_x = min(flatten_list(min_thks)) - interval_x
        max_y = interval_y - extreme_max % 10 + extreme_max
        min_y = extreme_min - extreme_min % 10
    else:
        [[max_x, min_x], [max_y, min_y]] = xy_limits

    if not GSF_ticks:
        last_scale = round(max_x/min_x, 2)
        middle_scale = round((last_scale - 1.0)/2 + 1, 2)
        quad_scale = round((last_scale - 1.0)*3/4 + 1, 2)
        ticks_GSF = [middle_scale, quad_scale, last_scale]
    else:
        ticks_GSF = GSF_ticks

    ticks_x = []
    for tck in ticks_GSF:
        ticks_x.append(max_x/tck)

    # print(thicknesses)
    # print(solutions)
    # print(kmax)
    for i in range(kmax):
        thks = thicknesses[i]
        sols = solutions[i]

        xmin = array(thks[0])
        xmax = array(thks[1])

        fmin = 100.0*array(sols[0])
        fmax = -100.0*array(sols[1])
        n = len(xmin)
        m = len(xmax)

        if simplified is True:
            ax.plot(xmin, fmin, markers[i], ls='-', markersize=6, color=colors[i], label=legends[i])
            ax.plot(xmax, fmax, markers[i], ls='-', markersize=6, color=colors[i])

        else:
            ax.plot(xmin, fmin, markers[i], ls='-', markersize=6, color='blue', label='minimum thrust '+legends[i])
            ax.plot(xmax, fmax, markers[i], ls='-', markersize=6, color='red', label='maximum thrust '+legends[i])

        if limit_state:
            x_, y_, _ = intersection_line_line_xy([[xmax[m-1], fmax[m-1]], [xmax[m-2], fmax[m-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
            extrapolation_maxy = [fmax[m-1], y_]
            extrapolation_miny = [fmin[n-1], y_]
            extrapolation_maxx = [xmax[m-1], x_]
            extrapolation_minx = [xmin[n-1], x_]
            ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black')
            ax.plot(extrapolation_maxx, extrapolation_maxy, '', ls='--', color=colors[i])
            ax.plot(extrapolation_minx, extrapolation_miny, '', ls='--', color=colors[i])
            ax.annotate(str(round(max_x/x_, 1)), (x_, y_), textcoords="offset points", xytext=(20, -5), ha='center', size=12)

        if fill:
            xmin_ = append(xmin, x_)
            xmax_ = append(xmax, x_)
            fmin_ = append(fmin, y_)
            fmax_ = append(fmax, y_)

            ax.fill(append(xmin_, xmax_[::-1]), append(fmin_, fmax_[::-1]), color=colors[i], alpha=0.2)
            # ax.fill_between(xmin, fmin, fmax, color=colors[i], alpha=0.2)
            # ax.fill_between(extrapolation_minx, extrapolation_miny, extrapolation_maxy, color=colors[i], alpha=0.2)

    ax2 = ax.twiny()
    print([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    print(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF])

    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
    ax.set_xlabel(r'thickness-over-span ($t/s$)', size=size_axis_label, labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)
    ax.set_ylabel(r'thrust/weight ($T_i/W_i$) [%]', size=size_axis_label, labelpad=8)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax.tick_params(axis='both', which='major', labelsize=size_axis_data)

    # plt.axvline(x=(x_ - max_x)/(min_x - max_x)*100, ls='--', color='black') # add    ymin=(max_y-y_)/(max_y-min_y)

    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    if show_legend is True:
        ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)
        box = ax.get_position()
        ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])

    if save:
        plt.savefig(save)
        print('Diagram saved at:', save)

    return plt

def box_limits(box):
    return [box.x0, box.y0, box.width*0.8, box.height]

size_plots = (10, 5)

size_axis_label = 12
size_axis_data = 12
size_legend = 12

# Plot for graph of curves and discretisations
series = []
x = []
y = []
styles = []
colours = []

radius = 5.0
bench = 0.042
discr = [4, 8, 12, 16, 20, 24]

# meridians = 12
series.append(r'$n_\mathrm{M} = 12$')
y.append(np.array([0.06670, 0.17474, 0.19103, 0.19473, 0.20282, 0.20311])/radius)

# meridians = 16
series.append(r'$n_\mathrm{M} = 16$')
y.append(np.array([0.06922, 0.17687, 0.19266, 0.19669, 0.20455, 0.20464])/radius)

# meridians = 20
series.append(r'$n_\mathrm{M} = 20$')
y.append(np.array([0.07037, 0.17786, 0.19341, 0.19759, 0.20534, 0.20534])/radius)

# meridians = 24
series.append(r'$n_\mathrm{M} = 24$')
y.append(np.array([0.07100, 0.17839, 0.19382, 0.19808, 0.20577, 0.20572])/radius)

# theoretical
series.append('benchmark')
bench_x = [4, 24]
y.append([bench, bench])

# for i in range(len(series)):
#     plt.plot(x[i], y[i], 'o-', label=series[i])


fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
ax = plt.axes()
ax.plot(discr, y[0], 'o-', label=series[0])
ax.plot(discr, y[1], 'o-', label=series[1])
ax.plot(discr, y[2], 'o-', label=series[2])
ax.plot(discr, y[3], 'o-', label=series[3])
ax.plot(bench_x, y[4], color='black', linestyle='dashed', label=series[4])
ax.legend(fontsize=size_legend)
ax.set_xlabel(r'number of parallels $(n_\mathrm{P})$', size=size_axis_label)
ax.set_ylabel(r'thickness-over-radius $(t/r)$', size=size_axis_label)
ax.annotate(r'benchmark $(t/r) = 0.042$', (sum(bench_x)/2, 0.042 + 0.001), textcoords="offset points", xytext=(sum(bench_x)/2, 0.042), ha='center')  # size =
ax.tick_params(labelsize=size_axis_data)
ax.set_xlim(4-1, 24+1)
ax.set_ylim(0, 0.05)
ax.set_xticks([4, 8, 12, 16, 20, 24])

box = ax.get_position()
ax.set_position(box_limits(box))

plt.show()

# ------------------- Plot of Dome

type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [20, 16]

folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

csv_file = os.path.join(folder, title + '_data.csv')
print(csv_file)
thicknesses, solutions = open_csv_row(csv_file, cut_last=False)
thickness_over_radius = [[],[]]
for el in thicknesses[0]:
    thickness_over_radius[0].append(el/radius)
    thickness_over_radius[1].append(el/radius)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thickness_over_radius, solutions, save=img_graph, fill=True, limit_state=False, x_label=r'thickness-over-radius ($t/r$)', GSF_ticks=[1.5, 2.0, 2.5, 3.0]).show()

# # ------------------- Plot of Amiens solution

# type_structure = 'crossvault'
# type_formdiagram = 'fan_fd'
# discretisation = 14
# file_name = 'amiens_internet'

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'max_n', file_name, type_structure, type_formdiagram, 'min_max')
# os.makedirs(folder, exist_ok=True)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_offset-method'

# csv_file = os.path.join(folder, title + '_data.csv')
# thicknesses, solutions = open_csv_row(csv_file, cut_last=False)
# print(thicknesses)

# img_graph = os.path.join(folder, title + '_diagram.pdf')
# img_graph = False
# diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False, GSF_ticks=[1.2, 1.4, 1.6]).show()

# ------------------- Plot of Crossvault graph

# path = '/Users/mricardo/compas_dev/me/min_thk/crossvault/study_crossvault_minthk.csv'

# discr = []
# cross_fd = [[], [], []]
# fan_fd = [[], [], []]
# radius = 10.0  #actually span
# # series = ['Analytical min thickness', 'Offset min thickness', 'Correction on normals']
# series = ['Orthogonal Form Diagram', 'Fan-like Form Diagram']

# with open(path) as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line = 0
#     for row in csv_reader:
#         if line > 0:
#             discr.append(int(row[0]))
#             for i in range(3):
#                 cross_fd[i].append(float(row[1 + 2*i])/radius)
#                 fan_fd[i].append(float(row[2 + 2*i])/radius)
#         line += 1

# print(cross_fd)
# print(fan_fd)
# lines = []
# fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
# lines += ax.plot(discr, cross_fd[0], 'o-', label=series[0])#, color='royalblue')
# # lines += ax.plot(discr, cross_fd[1], 's--', label=series[1])#, color='royalblue')
# # lines += ax.plot(discr, cross_fd[2], 'x--', label=series[2])#, color='royalblue')

# lines += ax.plot(discr, fan_fd[0], 'o-', label=series[1])#, color='mediumpurple')
# # lines += ax.plot(discr, fan_fd[1], 's--', label=series[1])#, color='mediumpurple')
# # lines += ax.plot(discr, fan_fd[2], 'x--', label=series[2])#, color='mediumpurple')

# # plt.plot(bench_x, y[4], color='black', linestyle='dashed', label=series[4])
# # ax.legend(lines[:2], ['line A', 'line B'],
# #           loc='upper right', frameon=False)


# from matplotlib.legend import Legend

# # ax = plt.axes()

# ax.legend(fontsize=size_legend)
# ax.tick_params(labelsize=size_axis_data)

# # ax.legend(lines[:3], series, loc='upper right', fontsize=size_legend, title='Orthogonal Diagram:', bbox_to_anchor=(1.35, 0.8))

# # leg = Legend(ax, lines[3:], series, frameon=False, loc='lower right', fontsize=size_legend, title='Fan-like Diagram:', bbox_to_anchor=(1.35, 0.25))
# # ax.add_artist(leg)

# box = ax.get_position()
# ax.set_position(box_limits(box))

# # leg = Legend(ax, lines[2:], ['line C', 'line D'],
# #              loc='lower right', frameon=False)

# ax.set_xlabel(r'discretisation $(n)$', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness-over-span $(t/s)$', size=size_axis_label, labelpad=8)
# # ax.annotate(r'benchmark $(t/R) = 0.042$', (sum(bench_x)/2, 0.042 + 0.001), textcoords="offset points", xytext=(sum(bench_x)/2, 0.042), ha='center')  # size =
# ax.set_xlim(discr[0]-1, discr[-1]+1)
# ax.set_ylim(0, 0.05)
# ax.set_xticks(discr)
# plt.show()

# # --------------------------   Plot spring angle interms of deg (discretisation n=14)

# # can plot in terms of A reading the file study_crossvault_A
# path = '/Users/mricardo/compas_dev/me/min_thk/crossvault/study_crossvault_deg_D=14_minthk.csv'

# A = []
# cross_fd = [[], [], [], []]
# fan_fd = [[], [], [], []]
# radius = 10.0
# series = ['analytical', 'strategy A', 'strategy B']

# with open(path) as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line = 0
#     for row in csv_reader:
#         if line > 0:
#             A.append(float(row[0]))
#             for i in range(4):
#                 cross_fd[i].append(float(row[1 + 2*i])/radius)
#                 fan_fd[i].append(float(row[2 + 2*i])/radius)
#         line += 1

# lines = []
# fig = plt.figure(figsize=size_plots)
# ax = plt.subplot(111)
# lines += ax.plot(A, cross_fd[0], 'o-', label='Orthogonal Form Diagram', color='C0')
# lines += ax.plot(A, fan_fd[0], 'o-', label='Fan-like Form Diagram', color='C1')

# ax.legend(fontsize=size_legend)
# ax.tick_params(labelsize=size_axis_data)
# box = ax.get_position()
# ax.set_position(box_limits(box))

# ax.set_xlabel(r'springing angle ($\beta$) [$^{\circ}$]', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness-over-span $(t/s)$', size=size_axis_label, labelpad=8)
# ax.set_xlim(A[0]-1, A[-1]+1)
# ax.set_ylim(0, 0.05)
# ax.set_xticks(A)
# plt.show()

# lines = []
# fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
# lines += ax.plot(A, cross_fd[0], 'o-', label=series[0], color='C0')
# lines += ax.plot(A, cross_fd[1], 's--', label=series[1], color='C2')
# lines += ax.plot(A, cross_fd[2], 'x--', label=series[2], color='C3')
# # lines += ax.plot(A[-6:], cross_fd[3][-6:], '.--', label=series[3], color='C4')

# ax.legend(fontsize=size_legend)
# ax.tick_params(labelsize=size_axis_data)
# box = ax.get_position()
# ax.set_position(box_limits(box))

# ax.set_xlabel(r'springing angle ($\beta$) [$^{\circ}$]', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness-over-span $(t/s)$', size=size_axis_label, labelpad=8)
# ax.set_xlim(A[0]-1, A[-1]+1)
# ax.set_ylim(0, 0.10)
# ax.set_xticks(A)
# plt.show()


# lines = []
# fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
# lines += ax.plot(A, fan_fd[0], 'o-', label=series[0], color='C1')
# lines += ax.plot(A, fan_fd[1], 's--', label=series[1], color='C2')
# lines += ax.plot(A, fan_fd[2], 'x--', label=series[2], color='C3')
# # lines += ax.plot(A[-6:], fan_fd[3][-6:], '.--', label=series[3], color='C4')

# ax.legend(fontsize=size_legend)
# ax.tick_params(labelsize=size_axis_data)
# box = ax.get_position()
# ax.set_position(box_limits(box))

# ax.set_xlabel(r'springing angle ($\beta$) [$^{\circ}$]', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness-over-span $(t/s)$', size=size_axis_label, labelpad=8)
# ax.set_xlim(A[0]-1, A[-1]+1)
# ax.set_ylim(0, 0.10)
# ax.set_xticks(A)
# plt.show()

# ------------- plot of the diagrams of thrust

xy_limits = [[0.050, 0.001], [210, 70]]
GSF_ticks = [2.0, 3.0, 4.0, 5.0]
legends = {'cross_fd': [r'orthogonal | $\beta=0^{\circ}$', r'orthogonal | $\beta=20^{\circ}$', r'orthogonal | $\beta=40^{\circ}$'], 'fan_fd': [r'fan-like | $\beta=0^{\circ}$', r'fan-like | $\beta=20^{\circ}$', r'fan-like | $\beta=40^{\circ}$']}
colors = {'cross_fd': ['C0', 'C0', 'C0'], 'fan_fd': ['C1', 'C1', 'C1']}  # These are C0 and C1 in HEX.
colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

type_structure = 'crossvault'
discretisation = 14
span = 10.0
for type_formdiagram in ['cross_fd', 'fan_fd']:
    thicknesses_all = []
    solutions_all = []
    for deg in [0, 20, 40]:
        folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'deg='+str(deg))
        os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_deg=' + str(deg)
        forms_address = os.path.join(folder, title)
        csv_file = os.path.join(folder, title + '_data.csv')
        thicknesses, solutions = open_csv_row(csv_file, cut_last=False)
        img_graph = None
        # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
        thicknesses_all.append(thicknesses)
        solutions_all.append(solutions)
        print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])
    folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
    title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    img_graph = os.path.join(folder_main, title_main + '_diagram.pdf')
    print(thicknesses_all)
    thickness_over_span = [[[], []], [[], []], [[], []]]
    for i in range(len(thicknesses_all)):
        for el in thicknesses_all[i][0]:
            print(el)
            thickness_over_span[i][0].append(el/span)
            thickness_over_span[i][1].append(el/span)

    img_graph = None
    diagram_of_multiple_thrust(thickness_over_span, solutions_all, legends[type_formdiagram], save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()
    # diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()
