from compas_plotters import MeshPlotter
from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from numpy import arange
from numpy import array
from numpy import append

from compas.geometry import intersection_line_line_xy

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'diagram_of_thrust',
    'diagram_of_multiple_thrust',
    'diagram_of_thrust_load_mult',
]


def diagram_of_thrust(dimension, min_sol, max_sol, limit_state=True, save=False):
    """ Plot a diagram of Thrusts based on the collected data from (n) points.

    Parameters
    ----------
    dimensions : list (n)
        Points with discretised solutions.
    min_sol : list (n)
        Adimensional thrust over weight for minimum thrust.
    max_sol : list (n)
        Adimensional thrust over weight for minimum thrust.

    Returns
    -------
    obj
        Plotter object.

    """
    xmin = xmax = array(dimension)
    fmin = 100.0 * array(min_sol)
    fmax = -100.0 * array(max_sol)
    n = len(xmin)
    x_, y_, _ = intersection_line_line_xy([[xmax[n-1], fmax[n-1]], [xmax[n-2], fmax[n-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
    print(x_, y_)
    xmax_ = append(xmax, x_)
    fmin_ = append(fmin, y_)
    fmax_ = append(fmax, y_)

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    interval_x = abs(dimension[0] - dimension[1])
    interval_y = 10
    max_x = max(xmax_)
    min_x = min(xmax_) - interval_x
    max_y = interval_y - max(fmax_) % 10 + max(fmax_)
    min_y = min(fmin_) - min(fmin_) % 10
    last_scale = round(max_x/min_x, 2)
    middle_scale = round((last_scale - 1.0)/2 + 1, 1)
    middle_x = max_x/middle_scale

    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='blue', label='minimum thrust')
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='red', label='maximum thrust')

    if limit_state:
        extrapolation_max = [fmax[n-1], y_]
        extrapolation_min = [fmin[n-1], y_]
        extrapolation_x = [xmax[n-1], x_]
        ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')
        ax.plot(extrapolation_x, extrapolation_max, '', ls='--', color='red')
        ax.plot(extrapolation_x, extrapolation_min, '', ls='--', color='blue')
        ax.annotate(str(round(x_, 2)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([0, 100*(max_x - middle_x)/(max_x-min_x), 100*(max_x - min(xmax))/(max_x-min_x), 100])
    ax2.set_xticklabels(['1.0', str(middle_scale), str(round(max_x/min(xmin), 2)), str(last_scale)], size=size_axis_data)
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

    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def diagram_of_multiple_thrust(dimensions, min_sols, max_sols, legends):
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

    kmax = len(dimensions)
    markers = ['o', '^', 's', 'D', 'x', '1', '2', '3']

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    fig = plt.figure(figsize=[12, 4])
    ax = fig.add_subplot(1, 1, 1)
    def flatten_list(l): return [item for sublist in l for item in sublist]
    extreme_max = -100*min(flatten_list(max_sols))
    extreme_min = 100*min(flatten_list(min_sols))

    interval_x = abs(flatten_list(dimensions)[0] - flatten_list(dimensions)[1])
    interval_y = 10
    max_x = max(flatten_list(dimensions))
    min_x = min(flatten_list(dimensions)) - interval_x
    max_y = interval_y - extreme_max % 10 + extreme_max
    min_y = extreme_min - extreme_min % 10
    last_scale = round(max_x/min_x, 2)
    middle_scale = round((last_scale - 1.0)/2 + 1, 1)
    middle_x = max_x/middle_scale

    print(max_x, min_x, max_y, min_y, min(min(dimensions)))

    for i in range(kmax):

        xmin = xmax = array(dimensions[i])
        fmin = 100.0*array(min_sols[i])
        fmax = -100.0*array(max_sols[i])
        print(xmin, fmin)
        print(xmax, fmax)

        # x_ = xmin[len(xmin)-1]
        # y_ = fmin[len(xmin)-1]

        ax.plot(xmin, fmin, markers[i], ls='-', markersize=6, color='blue', label='minimum thrust '+legends[i])
        ax.plot(xmax, fmax, markers[i], ls='-', markersize=6, color='red', label='maximum thrust '+legends[i])
        # ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([0, 100*(max_x - middle_x)/(max_x-min_x), 100])
    ax2.set_xticklabels(['1.0', str(middle_scale), str(last_scale)], size=size_axis_data)
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
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    return plt


def diagram_of_thrust_load_mult(dimension, min_sol, max_sol, limit_state=True, save=None):
    """ Plot a diagram of Thrusts to the problem of increasing a load multiplier based on the collected data from (n) points.

    Parameters
    ----------
    dimensions : list (n)
        Points with discretised solutions.
    min_sol : list (n)
        Adimensional thrust over weight for minimum thrust.
    max_sol : list (n)
        Adimensional thrust over weight for maximum thrust.
    limit_state : bool (true)
        If set to True, it will extrapolate data from the final 2 results and linearly extrapolate the limit state.
    save : str (optional)
        Path to save the plot, if wanted.

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

    ax1 = plt.axes()
    ax1.set_xlabel('lambda', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('T/W (%)', size=size_axis_label, weight='bold', labelpad=8)
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

    return pl


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    x = [0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14], [0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13]
    sol_min = [[0.32122703394577, 0.32864107398087034, 0.3363224906540916, 0.3442906794077901, 0.3525662450288712, 0.3611714918329792, 0.3701304796237221], [
        0.31802320817892216, 0.32546638697810754, 0.33317663927442764, 0.34117328843022754, 0.34947705652765554, 0.3581102296807935, 0.3670968859603074, 0.3764630724502445]]
    sol_max = [[-0.48442375197943505, -0.473126953779758, -0.4588129434920398, -0.44482604942259, -0.4312248478984081, -0.4154048423235026, -0.39113151140302316],
               [-0.49764114502249457, -0.4862719137722568, -0.47511606729242484, -0.4624620930999378, -0.448746719690402, -0.43539462435512494, -0.4223895081061913, -0.40971622703057375]]
    legends = ['a', 'b']
    diagram_of_multiple_thrust(x, sol_min, sol_max, legends).show()
