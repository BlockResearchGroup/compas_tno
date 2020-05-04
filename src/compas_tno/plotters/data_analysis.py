from compas_plotters import MeshPlotter
from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from numpy import arange
from numpy import array
from numpy import append
from numpy import linspace

from compas.geometry import intersection_line_line_xy

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import compas_tno
import csv
import os

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'diagram_of_thrust',
    'diagram_of_multiple_thrust',
    'diagram_of_thrust_load_mult',
    'surface_GSF_load_mult',
    'save_csv',
    'open_csv',
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
    limit_state : bool
        If yes, it interpolates the limit_state based on the two last solutions.
    limit_state : bool
        If yes, it saves the figure of the diagram.

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
    middle_scale = round((last_scale - 1.0)/2 + 1, 2)
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
        ax.annotate(str(round(x_, 3)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([0, 100*(max_x - middle_x)/(max_x-min_x), 100*(max_x - min(xmax_))/(max_x - min_x), 100])
    ax2.set_xticklabels(['1.0', str(middle_scale), str(round(max_x/min(xmax_), 2)), str(last_scale)], size=size_axis_data)
    ax1.set_xlabel('t/R', size=size_axis_label, weight='bold', labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('T/W [%]', size=size_axis_label, weight='bold', labelpad=8)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def diagram_of_multiple_thrust(dimensions, min_sols, max_sols, legends, simplified=True, limit_state=True, colors=None, xy_limits=None, save=None):
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

    if colors == None:
        colormap = plt.cm.coolwarm  # gist_ncar nipy_spectral, Set1, Paired coolwarm
        colors = [colormap(i) for i in linspace(0, 1.0, kmax)]

    size_axis_label = 14
    size_axis_data = 12
    size_legend = 14

    fig = plt.figure(figsize=[12, 4])
    ax = fig.add_subplot(1, 1, 1)
    def flatten_list(l): return [item for sublist in l for item in sublist]
    extreme_max = -100*min(flatten_list(max_sols))
    extreme_min = 100*min(flatten_list(min_sols))

    interval_x = abs(flatten_list(dimensions)[0] - flatten_list(dimensions)[1])
    interval_y = 20
    if xy_limits == None:
        max_x = max(flatten_list(dimensions))
        min_x = min(flatten_list(dimensions)) - interval_x
        max_y = interval_y - extreme_max % 10 + extreme_max
        min_y = extreme_min - extreme_min % 10
    else:
        [[max_x, min_x],[max_y, min_y]] = xy_limits

    last_scale = round(max_x/min_x, 2)
    middle_scale = round((last_scale - 1.0)/2 + 1, 1)
    middle_x = max_x/middle_scale

    # print(max_x, min_x, max_y, min_y, min(min(dimensions)))

    for i in range(kmax):

        xmin = xmax = array(dimensions[i])
        fmin = 100.0*array(min_sols[i])
        fmax = -100.0*array(max_sols[i])
        n = len(xmin)

        if simplified == True:
            ax.plot(xmin, fmin, markers[i], ls='-', markersize=6, color=colors[i], label=legends[i])
            ax.plot(xmax, fmax, markers[i], ls='-', markersize=6, color=colors[i])

        else:
            ax.plot(xmin, fmin, markers[i], ls='-', markersize=6, color='blue', label='minimum thrust '+legends[i])
            ax.plot(xmax, fmax, markers[i], ls='-', markersize=6, color='red', label='maximum thrust '+legends[i])

        if limit_state:
            x_, y_, _ = intersection_line_line_xy([[xmax[n-1], fmax[n-1]], [xmax[n-2], fmax[n-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
            extrapolation_max = [fmax[n-1], y_]
            extrapolation_min = [fmin[n-1], y_]
            extrapolation_x = [xmax[n-1], x_]
            ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black')
            ax.plot(extrapolation_x, extrapolation_max, '', ls='--', color=colors[i])
            ax.plot(extrapolation_x, extrapolation_min, '', ls='--', color=colors[i])
            ax.annotate(str(round(max_x/x_, 2)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    ax1 = plt.axes()
    ax2 = ax1.twiny()
    ax2.set_xticks([0, 100*(max_x - middle_x)/(max_x-min_x), 100])
    ax2.set_xticklabels(['1.0', str(middle_scale), str(last_scale)], size=size_axis_data)
    ax1.set_xlabel('t/R', size=size_axis_label, weight='bold', labelpad=8)
    ax2.set_xlabel('GSF', size=size_axis_label, weight='bold', labelpad=8)
    ax1.set_ylabel('T/W [%]', size=size_axis_label, weight='bold', labelpad=8)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax1.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    # ax1.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax1.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax1.tick_params(axis='both', which='major', labelsize=size_axis_data)

    # plt.axvline(x=(x_ - max_x)/(min_x - max_x)*100, ls='--', color='black') # add    ymin=(max_y-y_)/(max_y-min_y)

    box = ax.get_position()
    ax.set_position([box.x0*0.6, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

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
    sizes : list of lists (m-n)
        Points with discretised solutions.
    mins : list of lists (m-n)
        Adimensional thrust over weight for minimum thrust.
    maxs : list of lists (m-n)
        Adimensional thrust over weight for minimum thrust.
    legends : list (m)
        Legend of the (m) problems.
    save : str
        Path to save the diagram.

    Returns
    -------
    obj
        Plotter object.

    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    def flatten_list(l): return [item for sublist in l for item in sublist]
    def flatten_list_inv(l): return [-1*item for sublist in l for item in sublist]

    m = [len(sublist) for sublist in sizes]
    xs = array(flatten_list(sizes)+flatten_list(sizes))
    ys = array(flatten_list(mins)+flatten_list_inv(maxs))
    z_list = flatten_list([[legends[i]]*m[i] for i in range(len(legends))])
    zs = array(z_list + z_list)

    ax.plot_trisurf(xs, ys, zs, cmap=cm.coolwarm)
    ax.set_xlabel('t/R')
    ax.set_ylabel('T/W')
    ax.set_zlabel('px')

    if save:
        plt.savefig(save)

    return plt


def save_csv(dimension, min_sol, max_sol, limit_state=True, path=None, title=None):

    xmin = xmax = array(dimension)
    fmin = 100.0 * array(min_sol)
    fmax = -100.0 * array(max_sol)
    n = len(xmin)
    print(xmin)
    print(xmin.shape)
    if limit_state:
        x_, y_, _ = intersection_line_line_xy([[xmax[n-1], fmax[n-1]], [xmax[n-2], fmax[n-2]]], [[xmin[n-1], fmin[n-1]], [xmin[n-2], fmin[n-2]]])
        xmin = append(xmax, x_)
        fmin = append(fmin, y_)
        fmax = append(fmax, y_)

    if path is None:
        path = os.path.join(compas_tno.get('/csv/'),'test.csv')

    with open(path, mode='w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if title is not None:
            csv_writer.writerow([title])
        csv_writer.writerow(['X-Axis', 'min T/W', 'max T/W'])
        for i in range(len(fmin)):
            csv_writer.writerow([xmin[i], fmin[i], fmax[i]])

    return


def open_csv(path, cut_last=True):

    x = []
    min_thrust = []
    max_thrust = []
    line_total = 0

    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            line_total += 1

    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Title: {", ".join(row)}')
                line_count += 1
            elif line_count == 1:
                line_count += 1
                pass
            elif line_count < line_total - 1:
                x.append(float(row[0]))
                min_thrust.append(float(row[1])/100)
                max_thrust.append(-1 * float(row[2])/100)
                line_count += 1
        print(f'Processed {line_count} lines.')

    return x, min_thrust, max_thrust

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
    # diagram_of_multiple_thrust(x, sol_min, sol_max, legends).show()
    # save_csv(x[0], sol_min[0], sol_max[0], title='This is the example')


    x = [0.1, 0.096, 0.092, 0.088, 0.08399999999999999, 0.08, 0.076, 0.072, 0.068, 0.064, 0.06, 0.05600000000000001]
    ymin = [0.6733805282946738, 0.6768695278700919, 0.6801194890776822, 0.6830884474587336, 0.6857270723059036, 0.6879770995330894, 0.6897693561544739, 0.6910260162448293, 0.6916650604812549, 0.691535501586044, 0.6904910816351549, 0.6887132739276414]
    ymax = [-0.9339822060582105, -0.9153600350687701, -0.8967801607082414, -0.8749268981352708, -0.8527838337729975, -0.8308005117230084, -0.8090235909790843, -0.7873861336721817, -0.7658154927226042, -0.7442314431375009, -0.7225437965300713, -0.7004737016559901]
    diagram_of_thrust(x, ymin, ymax).show()
