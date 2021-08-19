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

    ax1 = plt.axes()
    ax2 = ax1.twiny()

    # ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    # ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF], size=size_axis_data)
    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)

    if x_label:
        ax1.set_xlabel(x_label, size=size_axis_label, labelpad=8)
    else:
        ax1.set_xlabel(r'thickness ($t$)', size=size_axis_label, labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)
    ax1.set_ylabel(r'thrust-over-weight ($T_i/W_i$) [%]', size=size_axis_label, labelpad=8)
    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # Check if this is necessary
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)
    ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def diagram_of_multiple_thrust(thicknesses, solutions, legends, simplified=True, limit_state=True, colors=None, xy_limits=None, GSF_ticks=False, fill=False, show_legend=True, save=None, markers=None, addpt=None, legendaddpt=None):
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

    interval_x = 0.05  # abs(flatten_list(min_thks)[0] - flatten_list(min_thks)[1])
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

        print(i)
        print(xmin)
        print(xmax)

        mult = 100
        sign = -1
        if sols[0][0] > 5.0:
            mult = 1
        if sols[1][0] > 0.0:
            sign = +1

        fmin = mult*array(sols[0])
        fmax = mult*sign*array(sols[1])
        n = len(xmin)
        m = len(xmax)
        print(fmin)
        print(fmax)

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

    if addpt:
        print('addpt')
        for i, (xi, yi) in enumerate(addpt):
            ax.plot(xi, yi, 'X', ls=' ', label=legendaddpt[i])

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    print([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    print(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF])

    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
    ax1.set_xlabel(r'thickness-over-span ($t/s$)', size=size_axis_label, labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)
    ax1.set_ylabel(r'thrust/weight ($T_i/W_i$) [%]', size=size_axis_label, labelpad=8)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

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

# ------------CROSS DIAGRAM------------
# ------------CROSS DIAGRAM------------
# ------------- plot of the diagrams of thrust

# xy_limits = [[0.50, 0.01], [160, 40]]
# GSF_ticks = [2.0, 3.0, 4.0, 5.0]
# legends = {'cross_fd': [r'orthogonal | $c=0$', r'orthogonal | $c=0.1$', r'orthogonal | $c=0.25$', r'orthogonal | $c=0.50$']}
# colors = {'cross_fd': ['C0', 'C1', 'C2', 'C3']}  # These are C0 and C1 in HEX.
# # colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# # colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

# type_structure = 'crossvault'
# discretisation = 14
# span = 10.0

# thicknesses_cross = {0: [], 0.5: [], 0.25: [], 0.10: []}
# solutions_cross = {0: [], 0.5: [], 0.25: [], 0.10: []}

# thicknesses_cross[0] = [[0.5, 0.45, 0.4, 0.35, 0.312962187], [0.5, 0.45, 0.4, 0.35, 0.312962187]]
# thicknesses_cross[0.1] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.258], [0.5, 0.45, 0.4, 0.35, 0.3, 0.258]]
# thicknesses_cross[0.25] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.184], [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.184]]
# thicknesses_cross[0.50] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08, 0.068], [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08, 0.068]]

# solutions_cross[0] = [[74.77527759, 76.80560202, 78.93177332, 81.16145847, 83.02732332], [101.4369571, 96.63002935, 91.66880597, 86.9492763, 83.02732332]]
# solutions_cross[0.1] = [[69.06129278, 71.08742835, 73.2219834, 75.47448344, 77.85555839, 82.99613507],
#                         [106.3914928, 102.2550684, 97.49670923, 92.52823744, 87.77523638, 82.99613507]]
# solutions_cross[0.25] = [[61.8501194, 63.6646636, 65.57630256, 67.59356665, 69.89140462, 72.72468194, 75.82010424, 83.02422223],
#                          [112.4785537, 108.7157143, 105.0077231, 100.7879474, 95.71380702, 90.75967139, 85.85373385, 83.02422223]]
# solutions_cross[0.50] = [[57.7043236, 58.76389089, 59.90347996, 61.08585247, 62.31355437, 63.58933168, 65.04157446, 68.78910449, 75.36704346, 78.36179929,
#                           80.76587164], [122.3514873, 118.2583604, 114.2776598, 110.385518, 106.5543606, 102.7497656, 97.83152217, 92.45423153, 86.52378904, 83.07284839, 80.76587164]]

# for type_formdiagram in ['cross_fd']:
#     thicknesses_all = []
#     solutions_all = []
#     for deg in [0, 0.1, 0.25, 0.50]:
#         thicknesses = thicknesses_cross[deg]
#         solutions = solutions_cross[deg]
#         img_graph = None
#         # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
#         thicknesses_all.append(thicknesses)
#         solutions_all.append(solutions)
#         print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])
#     folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
#     title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
#     img_graph = os.path.join(folder_main, title_main + '_diagram.pdf')
#     img_graph = None
#     print(thicknesses_all)
#     print(solutions_all)
#     # thickness_over_span = [[[], []], [[], []], [[], []]]
#     # for i in range(len(thicknesses_all)):
#     #     for el in thicknesses_all[i][0]:
#     #         print(el)
#     #         thickness_over_span[i][0].append(el/span)
#     #         thickness_over_span[i][1].append(el/span)

#     for series in solutions_all:
#         print(len(series[0]), len(series[1]))
#         for i in range(len(series[0])):
#             series[0][i] = series[0][i]/100
#             series[1][i] = series[1][i]/100
#     for series in thicknesses_all:
#         print(len(series[0]), len(series[1]))

#     print(thicknesses_all)
#     print(solutions_all)

#     diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph,
#                                fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()
#     # diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()



# ------------FAN DIAGRAM------------
# ------------FAN DIAGRAM------------
# ------------- plot of the diagrams of thrust

# xy_limits = [[0.50, 0.01], [160, 40]]
# GSF_ticks = [2.0, 3.0, 4.0, 5.0]
# legends = {'fan_fd': [r'fan-like | $c=0$', r'fan-like | $c=0.1$', r'fan-like | $c=0.25$', r'fan-like | $c=0.50$']}
# colors = {'fan_fd': ['C0', 'C1', 'C2', 'C3']}  # These are C0 and C1 in HEX.
# # colors = {'fan_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# # colors = {'fan_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

# type_structure = 'crossvault'
# discretisation = 14
# span = 10.0

# thicknesses_fan = {0: [], 0.5: [], 0.25: [], 0.10: []}
# solutions_fan = {0: [], 0.5: [], 0.25: [], 0.10: []}

# thicknesses_fan[0] = [[0.5, 0.462], [0.5, 0.462]]
# thicknesses_fan[0.1] = [[0.5, 0.45, 0.4, 0.381], [0.5, 0.45, 0.4, 0.381]]
# thicknesses_fan[0.25] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.284], [0.5, 0.45, 0.4, 0.35, 0.3, 0.284]]
# thicknesses_fan[0.50] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.133], [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.133]]

# solutions_fan[0] = [[85.39414706, 87.84789785], [91.58440479, 87.84789785]]
# solutions_fan[0.1] = [[77.38786732, 79.82182979, 82.92104332, 85.87888901], [98.61784773, 96.65808853, 91.82967859, 85.87888901]]
# solutions_fan[0.25] = [[69.30562706, 71.5898352, 74.01799122, 76.63720002, 79.43758131, 81.79817597], [108.6253247, 103.3151941, 98.18929977, 93.89348244, 88.71289272, 81.79817597]]
# solutions_fan[0.50] = [[58.32042403, 60.26863607, 62.34047242, 64.54856747, 65.92274897, 68.40896737, 72.20575127, 79.29859657], [121.5, 114.673246, 110.7220752, 105.7359792, 100.3871226, 95.27980421, 90.68620846, 79.29859657]]

# for type_formdiagram in ['fan_fd']:
#     thicknesses_all = []
#     solutions_all = []
#     for deg in [0, 0.1, 0.25, 0.5]:
#         thicknesses = thicknesses_fan[deg]
#         solutions = solutions_fan[deg]
#         img_graph = None
#         # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
#         thicknesses_all.append(thicknesses)
#         solutions_all.append(solutions)
#         print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])
#     folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
#     title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
#     img_graph = os.path.join(folder_main, title_main + '_diagram.pdf')
#     img_graph = None
#     print(thicknesses_all)
#     print(solutions_all)
#     # thickness_over_span = [[[], []], [[], []], [[], []]]
#     # for i in range(len(thicknesses_all)):
#     #     for el in thicknesses_all[i][0]:
#     #         print(el)
#     #         thickness_over_span[i][0].append(el/span)
#     #         thickness_over_span[i][1].append(el/span)

#     for series in solutions_all:
#         print(len(series[0]), len(series[1]))
#         for i in range(len(series[0])):
#             series[0][i] = series[0][i]/100
#             series[1][i] = series[1][i]/100
#     for series in thicknesses_all:
#         print(len(series[0]), len(series[1]))

#     print(thicknesses_all)
#     print(solutions_all)

#     diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph,
#                                fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()
#     # diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()


# ------------------------------------------
# ------------ CROSS-FD ------------
# ------------- plot of the diagrams of thrust DOUBLE CURVATURE

xy_limits = [[0.50, 0.01], [160, 40]]
GSF_ticks = [2.0, 3.0, 4.0, 5.0]
legends = {'fan_fd': [r'fan-like | $c=0$', r'fan-like | $c=0.1$'], 'cross_fd': [r'orthogonal | $c=0$', r'orthogonal | $c=0.1$']}
colors = {'fan_fd': ['C0', 'C1'], 'cross_fd': ['C2', 'C3']}  # These are C0 and C1 in HEX.
legend_all = []
# colors = {'fan_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# colors = {'fan_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

type_structure = 'crossvault'
discretisation = 14
span = 10.0

thicknesses = {'cross_fd': {0: [], 0.10: []}, 'fan_fd': {0: [], 0.10: []}}
solutions = {'cross_fd': {0: [], 0.10: []}, 'fan_fd': {0: [], 0.10: []}}

thicknesses['fan_fd'][0] = [[0.5, 0.469], [0.5, 0.469]]
thicknesses['fan_fd'][0.1] = [[0.5, 0.45, 0.4, 0.367], [0.5, 0.45, 0.4, 0.367]]

thicknesses['cross_fd'][0] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.285], [0.5, 0.45, 0.4, 0.35, 0.3, 0.285]]
thicknesses['cross_fd'][0.1] = [[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.218], [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.218]]


solutions['fan_fd'][0] = [[85.39414706, 87.84789785], [91.58440479, 87.84789785]]
solutions['fan_fd'][0.1] = [[77.38786732, 79.82182979, 82.92104332, 85.87888901], [98.61784773, 96.65808853, 91.82967859, 85.87888901]]

solutions['cross_fd'][0] = [[64.34443, 65.93296, 67.58872, 69.32432, 71.14970, 73.17908], [86.061436, 83.551783, 81.000579, 78.263841, 75.350514, 73.179083]]
solutions['cross_fd'][0.1] = [[58.92723779, 60.40426829, 61.9776272, 63.63727805, 65.40329015, 67.31489417, 71.79642535], [89.4634744, 86.9090431, 84.34693797, 81.7439412, 79.10048825, 76.07661758, 71.79642535]]

thicknesses_all = []
solutions_all = []
colors_all = []
# print(thicknesses)
for type_formdiagram in ['fan_fd', 'cross_fd']:
    i = 1
    for deg in [0.1]:
        thicknesses_i = thicknesses[type_formdiagram][deg]
        solutions_i = solutions[type_formdiagram][deg]
        color = colors[type_formdiagram]
        legend_all.append(legends[type_formdiagram][i])
        img_graph = None
        # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
        thicknesses_all.append(thicknesses_i)
        solutions_all.append(solutions_i)
        colors_all.append(color[i])
        print(type_formdiagram, deg, -solutions_i[1][0]/solutions_i[0][0])
        i += 1
    folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
    title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    img_graph = os.path.join(folder_main, title_main + '_diagram.pdf')
    img_graph = None
    print(thicknesses_all)
    print(solutions_all)
    # thickness_over_span = [[[], []], [[], []], [[], []]]
    # for i in range(len(thicknesses_all)):
    #     for el in thicknesses_all[i][0]:
    #         print(el)
    #         thickness_over_span[i][0].append(el/span)
    #         thickness_over_span[i][1].append(el/span)

    # for series in solutions_all:
    #     print(len(series[0]), len(series[1]))
    #     for i in range(len(series[0])):
    #         series[0][i] = series[0][i]/100
    #         series[1][i] = series[1][i]/100
    # for series in thicknesses_all:
    #     print(len(series[0]), len(series[1]))

    x_y_add = [[0.372, 75.0], [0.296, 70.6], [0.325, 73.9], [0.267, 69.8]]

    legendaddpt = ['B1', 'B2', 'D3', 'D4']

    diagram_of_multiple_thrust(thicknesses_all, solutions_all, legend_all, save=img_graph,
                               fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors_all, addpt=x_y_add, legendaddpt=legendaddpt).show()
    # diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends[type_formdiagram], save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors[type_formdiagram]).show()

