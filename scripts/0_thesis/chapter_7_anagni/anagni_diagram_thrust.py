import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from compas_tno.utilities import open_csv_row
# from compas_tno.plotters import diagram_of_thrust
# from compas_tno.utilities import diagram_of_multiple_thrust
from matplotlib.ticker import FormatStrFormatter
from numpy import arange
from numpy import array
from numpy import append
from numpy import linspace
from compas.geometry import intersection_line_line_xy
import csv


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
    plt.rcParams["font.family"] = "Arial"
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
    ax.set_xlabel(r'thickness [m]', size=size_axis_label, labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)
    ax.set_ylabel(r'$T/W$ [%]', size=size_axis_label, labelpad=8)
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


size_axis_label = 12
size_axis_data = 10
size_legend = 10

xy_limits = [[0.25, 0.0], [110, 50]]
GSF_ticks = [1.5, 2.0, 3.0, 5.0]
deg = 20
Rs = [6.1147, 7.08, 7.9422]
# legends = {'cross_fd': [r'orthogonal | $\beta=20$'], 'fan_fd': [r'fan-like | $\beta=20$'], 'topology-crossbraced': [r'braced | $\beta=20$']}
# colors = {'cross_fd': ['C0'], 'fan_fd': ['C1'], 'topology-crossbraced': ['C2']}  # These are C0 and C1 in HEX.
colors = ['C0', 'C1', 'C2']
legends = ['diagram (c)', 'diagram (d)', 'diagram (e)']
# colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

type_structure = 'pointed_crossvault'
discretisation = 14
option = ['A', 'B', 'C']

thks = {}
mins = {}
maxs = {}

# data

thks['diagram (c)'] = [0.25, 0.23, 0.21, 0.19, 0.17, 0.15, 0.13, 0.11455718997617]
mins['diagram (c)'] = [0.688494419242597, 0.707439247664863, 0.727534177946957, 0.749080111336226, 0.772550677854392, 0.803126199066912, 0.836636093153821, 0.89078079260606]
maxs['diagram (c)'] = [-1.06881386721609, -1.0477352323341, -1.02676379258101, -1.0057044838431, -0.98485834282685, -0.965047576246748, -0.935991542194311, -0.89078079260606]

thks['diagram (d)'] = [0.25, 0.23, 0.21, 0.19, 0.169999999999999, 0.15, 0.13, 0.105562104633249]
mins['diagram (d)'] = [0.606328653561807, 0.619462163668277, 0.633562215501158, 0.649070836241779, 0.674595270789398, 0.704151302961082, 0.736515813454131, 0.821810138268751]
maxs['diagram (d)'] = [-1.05933211914047, -1.02877851473101, -1.00053321387393, -0.973136593263799, -0.945405898188295, -0.917681710399176, -0.888359831883978, -0.821810138268751]

thks['diagram (e)'] = [0.25, 0.23, 0.21, 0.19, 0.169999999999999, 0.15, 0.13, 0.109999999999999, 0.096795447487448]
mins['diagram (e)'] = [0.697577440050666, 0.716119852659212, 0.735699353282211, 0.756414904233467, 0.778364697202226, 0.801657564835727, 0.829304030987453, 0.863645296466506, 0.888729956011649]
maxs['diagram (e)'] = [-1.06397845463424, -1.04239020331161, -1.0215335966826, -1.00069774233659, -0.979711405761614, -0.95933260262061, -0.939480944895213, -0.917643485722287, -0.888729955930222]

i_ = 0

thicknesses_all = []
solutions_all = []

for legend in legends:

    # for type_formdiagram in ['fan_fd', 'cross_fd', 'topology-crossbraced']:
    #     folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(R), 'min_thk', 'deg='+str(deg))
    #     # os.makedirs(folder, exist_ok=True)
    #     title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    #     forms_address = os.path.join(folder, title)
    #     csv_file = os.path.join(folder, title + '_data.csv')
    #     thicknesses, solutions = open_csv_row(csv_file, cut_last=True)
    #     img_graph = None
    #     adim_thk = [[], []]
    #     for i in range(len(thicknesses)):
    #         for j in range(len(thicknesses[i])):
    #             adim_thk[i].append(thicknesses[i][j]/10.0)
    #     print(adim_thk)
    #     # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
    #     thicknesses_all.append(adim_thk)
    #     solutions_all.append(solutions)
    #     print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])

    thicknesses_all.append([thks[legend], thks[legend]])
    solutions_all.append([mins[legend], maxs[legend]])

    print('thicknesses', thicknesses_all)
    print('solutionss', solutions_all)

diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors).show()

