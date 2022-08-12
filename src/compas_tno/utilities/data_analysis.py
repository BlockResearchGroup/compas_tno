
from numpy import arange
from numpy import array
from numpy import append
from numpy import linspace

from compas.geometry import intersection_line_line_xy

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm

import compas_tno
import json
import csv
import os


__all__ = [
    'diagram_of_thrust',
    'diagram_of_multiple_thrust',
    'diagram_of_thrust_load_mult',
    'surface_GSF_load_mult',
    'save_csv_row',
    'open_csv_row',
    'interpolate_min_thk',
    'filter_min_thk',
    'lookup_folder',
    'save_pointcloud'
]


def diagram_of_thrust(thicknesses, solutions, limit_state=True, fill=False, xy_limits=None, GSF_ticks=False, save=False):
    """ Plot a diagram of Thrusts based on the collected data from (n) points.

    Parameters
    ----------
    thicknesses : list of lists [[n],[m]]
        Points with discretised solutions for minimum/maximum thrust.
    solutions : list of lists [[n],[m]]
        Adimensional thrust over weight for minimum/maximum thrust.
    limit_state : bool, optional
        If yes, it interpolates the limit_state based on the two last solutions.
        The default value is ``True``.
    fill : bool, optional
        If yes, it fills bewteen the curves of min and max.
        The default value is ``False``.
    xy_limits : list of lists, optional
        Freezes the axis of the graphs.
        The default value is ``None`` in which case these limits are calculated automarically
    GSF_ticks : bool, optional
        Add the ticks regarding the GSF to the top axis.
        The default value is ``False``.
    save : str, optional
        Path to save the graph.
        The default value is ``False``.

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
    fmax = 100.0 * array(max_sol)  # add minus signe
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
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='blue', label='minimum thrust')
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='red', label='maximum thrust')

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

    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF], size=size_axis_data)

    ax1.set_xlabel('thickness/radius', size=size_axis_label, weight='bold', labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('thrust/weight [%]', size=size_axis_label, weight='bold', labelpad=8)
    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # Check if this is necessary
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

    # box = ax.get_position()
    # ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)
    ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def diagram_of_multiple_thrust(thicknesses, solutions, legends, simplified=True, limit_state=True, colors=None, xy_limits=None, GSF_ticks=False,
                               fill=False, show_legend=True, save=None, markers=None):
    """ Plot a diagram of Thrusts based on the collected data on (m) problems each with (n) points.

    Parameters
    ----------
    thicknesses : list of lists [[n],[m]]
        Points with discretised solutions for minimum/maximum thrust.
    solutions : list of lists [[n],[m]]
        Adimensional thrust over weight for minimum/maximum thrust.
    legends : list
        A list with the name of each data entry
    simplified : bool, optional
        If yes, the maximum and minimum thrust series have the same color.
        The default value is ``True``.
    limit_state : bool, optional
        If yes, it interpolates the limit_state based on the two last solutions.
        The default value is ``True``.
    colors : list, optional
        If added, this define the HEX color for each of the diagrams to plot.
        The default value is ``True``.
    xy_limits : list of lists, optional
        Freezes the axis of the graphs.
        The default value is ``None`` in which case these limits are calculated automarically
    GSF_ticks : bool, optional
        Add the ticks regarding the GSF to the top axis.
        The default value is ``False``.
    fill : bool, optional
        If yes, it fills bewteen the curves of min and max.
        The default value is ``False``.
    show_legend : bool, optional
        Defines if the legend is plotted in the graph.
        The default value is ``True``.
    markers : list, optional
        If added, this define the markers for each diagram.
        The default value is ``True``.
    save : str, optional
        Path to save the graph.
        The default value is ``False``.

    Returns
    -------
    obj
        Plotter object.

    """

    kmax = len(legends)

    if markers is None:
        markers = ['o', '^', 's', 'D', 'x', '1', '2', '3', 'v', 'p', '4', '8']

    if colors is None:
        colormap = plt.cm.get_cmap('coolwarm')  # gist_ncar nipy_spectral, Set1, Paired coolwarm
        colors = [colormap(i) for i in linspace(0, 1.0, kmax)]

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    fig = plt.figure(figsize=[15, 5])
    ax = fig.add_subplot(1, 1, 1)

    interval_x = 0.05  # abs(flatten_list(min_thks)[0] - flatten_list(min_thks)[1])
    interval_y = 20

    if xy_limits is None:
        min_sols = solutions[0][0]  # this has to find the maximum from all the input series
        max_sols = solutions[0][1]
        min_thks = thicknesses[0][1]
        def flatten_list(li): return [item for sublist in li for item in sublist]
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

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    print([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    print(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF])

    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
    ax1.set_xlabel('thickness/span', size=size_axis_label, weight='bold', labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('thrust/weight [%]', size=size_axis_label, weight='bold', labelpad=8)
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

        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)
        # box = ax.get_position()
        # ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])

    if save:
        plt.savefig(save)
        print('Diagram saved at:', save)

    return plt


def diagram_of_thrust_load_mult(dimension, min_sol, max_sol, limit_state=True, fill=False, xy_limits=None, save=None):
    """ Plot a diagram of Thrusts to the problem of increasing a load multiplier based on the collected data from (n) points.

    Parameters
    ----------
    dimensions : list (n)
        Points with discretised solutions.
    min_sol : list (n)
        Adimensional thrust over weight for minimum thrust.
    max_sol : list (n)
        Adimensional thrust over weight for maximum thrust.
    limit_state : bool, optional
        If yes, it interpolates the limit_state based on the two last solutions.
        The default value is ``True``.
    fill : bool, optional
        If yes, it fills bewteen the curves of min and max.
        The default value is ``False``.
    xy_limits : list of lists, optional
        Freezes the axis of the graphs.
        The default value is ``None`` in which case these limits are calculated automarically.
    save : str, optional
        Path to save the graph.
        The default value is ``False``.

    Returns
    -------
    obj
        Plotter object.

    """

    xmin = xmax = array(dimension)
    fmin = 100.0 * array(min_sol)
    fmax = -100.0 * array(max_sol)
    n = len(xmin)

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    interval_x = abs(dimension[0] - dimension[1])
    interval_y = 10
    if xy_limits:
        [[max_x, min_x], [max_y, min_y]] = xy_limits
    else:
        max_x = max(xmax) + interval_x
        min_x = min(xmax)
        max_y = interval_y - max(fmax) % 10 + max(fmax)
        min_y = min(fmin) - min(fmin) % 10

    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='blue', label='minimum thrust')
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='red', label='maximum thrust')

    if limit_state:
        x_, y_, _ = intersection_line_line_xy([[xmax[n-1], fmax[n-1]], [xmax[n-2], fmax[n-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
        extrapolation_max = [fmax[n-1], y_]
        extrapolation_min = [fmin[n-1], y_]
        extrapolation_x = [xmax[n-1], x_]
        ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')
        ax.plot(extrapolation_x, extrapolation_max, '', ls='--', color='red')
        ax.plot(extrapolation_x, extrapolation_min, '', ls='--', color='blue')
        ax.annotate(str(round(x_, 1)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    if fill:
        ax.fill_between(xmin, fmin, fmax, color='grey', alpha=0.2)
        ax.fill_between(extrapolation_x, extrapolation_min, extrapolation_max, color='grey', alpha=0.2)

    ax1 = plt.axes()
    ax1.set_xlabel('lambda', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('T/W [%]', size=size_axis_label, weight='bold', labelpad=8)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(min_x, max_x)
    ax1.set_ylim(min_y, max_y)
    ax1.set_xticks(arange(min_x, max_x + interval_x, interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)
    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def surface_GSF_load_mult(sizes, mins, maxs, legends, save=None):
    """ Plot a diagram of Thrusts based on the collected data on (m) problems each with (n) points.

    Parameters
    ----------
    sizes : list of lists [m-n]
        Points with discretised solutions.
    mins : list of lists [m-n]
        Adimensional thrust over weight for minimum thrust.
    maxs : list of lists [m-n]
        Adimensional thrust over weight for minimum thrust.
    legends : list [m]
        Legend of the (m) problems.
    save : str, optional
        Path to save the graph.
        The default value is ``False``.

    Returns
    -------
    obj
        Plotter object.

    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    def flatten_list(li): return [item for sublist in li for item in sublist]
    def flatten_list_inv(li): return [-1*item for sublist in li for item in sublist]

    m = [len(sublist) for sublist in sizes]
    xs = array(flatten_list(sizes)+flatten_list(sizes))
    ys = array(flatten_list(mins)+flatten_list_inv(maxs))
    z_list = flatten_list([[legends[i]]*m[i] for i in range(len(legends))])
    zs = array(z_list + z_list)

    ax.plot_trisurf(xs, ys, zs, cmap=cm.get_cmap('coolwarm'))
    ax.set_xlabel('t/R')
    ax.set_ylabel('T/W')
    ax.set_zlabel('px')

    if save:
        plt.savefig(save)

    return plt


def save_csv_row(thicknesses, solutions, limit_state=True, path=None, title=None):
    """Save a CSV file from the routine with several min/max thrust computations

    Parameters
    ----------
    thicknesses : list
        The list of the thickness computed for min/max thrust
    solutions : [list, list]
        The list of solutions of min/max solutions
    limit_state : bool, optional
        If limit state is included, by default True
    path : str, optional
        Full path to save the file, by default None
    title : str, optional
        Title for the analysis, by default None

    Returns
    -------
    None
        File is saved
    """

    xmin = array(thicknesses[0])
    xmax = array(thicknesses[1])
    fmin = 100.0 * array(solutions[0])
    fmax = -100.0 * array(solutions[1])
    n = len(xmin)
    m = len(xmax)

    if limit_state:
        x_, y_, _ = intersection_line_line_xy([[xmax[m-1], fmax[m-1]], [xmax[m-2], fmax[m-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
        # print('Calculated Limit State: (x, y):', x_, y_)
        xmin = append(xmin, x_)
        xmax = append(xmax, x_)
        fmin = append(fmin, y_)
        fmax = append(fmax, y_)

    if path is None:
        path = os.path.join(compas_tno.get('/csv/'), 'test.csv')

    with open(path, mode='w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if title:
            csv_writer.writerow([title])
        csv_writer.writerow(['X-Min'] + list(xmin))
        csv_writer.writerow(['X-Max'] + list(xmax))
        csv_writer.writerow(['min T/W'] + list(fmin))
        csv_writer.writerow(['max T/W'] + list(fmax))

    return


def open_csv_row(path, cut_last=True, printout=True):
    """Open a CSV file from the routine with several min/max thrust computations

    Parameters
    ----------
    path : str, optional
        Full path to save the file, by default None
    cut_last : bool, optional
        If should cut last element, by default True
    printout : bool, optional
        if should print in the screen, by default True

    Returns
    -------
    thicknesses : list
        The list of the thickness computed for min/max thrust
    solutions : [list, list]
        The list of solutions of min/max solutions
    """

    xmin = []
    xmax = []
    min_thrust = []
    max_thrust = []

    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        if cut_last is False:
            for row in csv_reader:
                if line_count == 0:
                    if printout:
                        print(f'Title-notcut: {", ".join(row)}')
                elif line_count == 1:
                    xmin.append(row[1:])
                elif line_count == 2:
                    xmax.append(row[1:])
                elif line_count == 3:
                    min_thrust.append(row[1:])
                elif line_count == 4:
                    max_thrust.append(row[1:])
                line_count += 1
        else:
            for row in csv_reader:
                if line_count == 0:
                    if printout:
                        print(f'Title: {", ".join(row)}')
                elif line_count == 1:
                    cut_length = len(row[1:])
                    xmin.append(row[1:cut_length])
                elif line_count == 2:
                    cut_length = len(row[1:])
                    xmax.append(row[1:cut_length])
                elif line_count == 3:
                    cut_length = len(row[1:])
                    min_thrust.append(row[1:cut_length])
                elif line_count == 4:
                    cut_length = len(row[1:])
                    max_thrust.append(row[1:cut_length])
                line_count += 1

    xmin = [float(i) for i in xmin[0]]
    xmax = [float(i) for i in xmax[0]]
    min_thrust = [float(i)/100 for i in min_thrust[0]]
    max_thrust = [-1 * float(i)/100 for i in max_thrust[0]]

    return [xmin, xmax], [min_thrust, max_thrust]


def interpolate_min_thk(sizes, solutions):
    """Interpolate min thickness from sizes and solutions computed

    Parameters
    ----------
    sizes : list
        Sizes computed
    solutions : list
        SOlutions of min/max optimisation

    Returns
    -------
    x_
        Point to add to the Stability domain with the interpolated value
    """

    xmin = array(sizes[0])
    xmax = array(sizes[1])
    fmin = 100.0 * array(solutions[0])
    fmax = -100.0 * array(solutions[1])
    n = len(xmin)
    m = len(xmax)
    x_, _, _ = intersection_line_line_xy([[xmax[m-1], fmax[m-1]], [xmax[m-2], fmax[m-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])

    return x_


def lookup_folder(folder):
    """Lookup files in a folder

    Parameters
    ----------
    folder : str
        Path to the folder

    Returns
    -------
    files_dict
        Dictionary with the information about the files
    """

    files = os.listdir(folder)
    # print(files)
    jsons = []
    for f in files:
        extension = f.split('.')[-1]
        if extension == 'json':
            jsons.append('.'.join(f.split('.')[:-1]))

    files_dict = {}

    for title in jsons:

        title_split = title.split('_')
        thk = float(title_split[-1])
        type_structure = '_'.join(title_split[:2])
        type_formdiagram = title_split[2]

        if type_formdiagram == 'cross' or type_formdiagram == 'fan':
            type_formdiagram = type_formdiagram + '_fd'

        if 'max' in title_split:
            objective = 'max'
        elif 'min' in title_split:
            objective = 'min'
        elif 'lp' in title_split:
            objective = 'lp'
        elif 'bestfit' in title_split:
            objective = 'bestfit'
        else:
            objective = 'other'
        if 'sag' in title:
            sag = int(title.split('sag_')[-1].split('_')[0])
        else:
            sag = False
        if 'smooth' in title:
            smooth = True
        else:
            smooth = False

        discr_sag_smooth = title_split[title_split.index('discr')+1]

        if not sag and not smooth:
            discr = int(discr_sag_smooth)
        elif sag:
            discr = int(discr_sag_smooth.split('sag')[0])
        else:
            discr = int(discr_sag_smooth.split('smooth')[0])

        data_file = {
            'thk': thk,
            'discretisation': discr,
            'type_structure': type_structure,
            'type_formdiagram': type_formdiagram,
            'objective': objective,
            'sag': sag,
            'smooth': smooth,
        }

        files_dict[title] = data_file

    return files_dict


def filter_min_thk(files_dict, filters=None):
    """Filter structure with minimum thickness in a folder

    Parameters
    ----------
    files_dict : dict
        Dictionary of all files
    filters : dict, optional
        Filters applied, by default None

    Returns
    -------
    limit_form_min, limit_form_max
        The limit values
    """

    limit_thk_min = 100
    limit_thk_max = 100
    limit_form_min = None
    limit_form_max = None

    for key, values in files_dict.items():
        proceed = True
        if filters:
            for filter in filters:
                if filters[filter] == values[filter]:
                    pass
                else:
                    proceed = False
                    break
        if proceed:
            if values['objective'] == 'min':
                if values['thk'] < limit_thk_min:
                    limit_thk_min = values['thk']
                    limit_form_min = {key: values}
            if values['objective'] == 'max':
                if values['thk'] < limit_thk_max:
                    limit_thk_max = values['thk']
                    limit_form_max = {key: values}

    return limit_form_min, limit_form_max


def save_pointcloud(points_lb, points_ub, json_path):
    """Save pointcloud to a JSON file

    Parameters
    ----------
    points_lb : list
        List of points of the intrados
    points_ub : list
        List of points of the extrados
    json_path : str
        Path to save the file
    """

    data = {'UB': {}, 'LB': {}}

    i = 0
    for pt in points_lb:
        data['LB'][i] = pt
        i += 1

    print('Found {0} points in LB'.format(i))

    i = 0
    for pt in points_ub:
        data['UB'][i] = pt
        i += 1

    print('Found {0} points in UB'.format(i))

    with open(json_path, 'w') as outfile:
        json.dump(data, outfile)

    return


# def prune_data(sizes, solutions):

#     xmin = sizes[0]
#     xmax = sizes[1]
#     fmin = solutions[0]
#     fmax = solutions[1]
#     n = len(xmin)
#     m = len(xmax)

#     xmax_ = []
#     xmin_ = []
#     fmax_ = []
#     fmin_ = []

#     if n == m:
#         return sizes, solutions
#     elif n > m and xmin[-1] == xmax[-1]:
#         xmax_ = xmax
#         fmax_ = fmax
#         j = 0
#         for i in range(m):
#             while xmax[i] != xmin[j]:
#                 j = j+1
#             if xmax[i] == xmin[j]:
#                 xmin_.append(xmin[j])
#                 fmin_.append(fmin[j])
#                 j = j+1
#         return [xmin_, xmax_], [fmin_, fmax_]
#     elif n < m and xmin[-1] == xmax[-1]:
#         xmin_ = xmin
#         fmin_ = fmin
#         j = 0
#         for i in range(n):
#             while xmin[i] != xmax[j]:
#                 j = j+1
#             if xmin[i] == xmax[j]:
#                 xmax_.append(xmax[j])
#                 fmax_.append(fmax[j])
#         return [xmin_, xmax_], [fmin_, fmax_]
#     else:
#         print('Last optimisation not with same thickness. Can not prune')
#         print([xmin, xmax], [fmin, fmax])
#         if m > n:
#             diff = m - n
#             print(diff)
#             for i in range(diff):
#                 xmax.pop(m-3)
#                 fmax.pop(m-3)
#             print([xmin, xmax], [fmin, fmax])
#             return [xmin, xmax], [fmin, fmax]
#         else:
#             diff = n - m
#             print(diff)
#             for i in range(diff):
#                 xmin.pop(m-3)
#                 fmin.pop(m-3)
#             print([xmin, xmax], [fmin, fmax])
#             return [xmin, xmax], [fmin, fmax]
#     return sizes, solutions
