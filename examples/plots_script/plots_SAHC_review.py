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
                fmin.append(float(row[1]))
                fmax.append(float(row[3])*-1.0)
                line_count += 1
        print(f'Processed {line_count} lines.')


    return [xmin, xmax], [fmin, fmax]


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
    fmax = -100.0 * array(max_sol)  # decide or not to invert
    # fmax = 100.0 * array(max_sol)
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
    print(xmax)
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

    # ax1 = plt.axes()
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

def diagram_of_multiple_thrust(thicknesses, solutions, legends, simplified=True, limit_state=True, colors=None, xy_limits=None, GSF_ticks=False, fill=False, show_legend=True, x_label=None, save=None, markers=None):
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

    print(legends)
    kmax = len(legends)
    print('kmax aaaaaaa', kmax)

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

    print(thicknesses)
    print(solutions)
    print(kmax)
    print('kmax', kmax)
    for i in range(kmax):
        print(i)
        thks = thicknesses[i]
        sols = solutions[i]

        xmin = array(thks[0])
        xmax = array(thks[1])

        fmin = 100.0*array(sols[0])
        fmax = -100.0*array(sols[1])
        n = len(xmin)
        m = len(xmax)

        if simplified is True:
            print(i)
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

    # ax1 = plt.axes()
    ax2 = ax.twiny()
    print([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    print(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF])

    ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
    ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
    if x_label:
        ax.set_xlabel(x_label, size=size_axis_label, labelpad=8)
    else:
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
radius = 5.0


# ------------------- Plot of Dome Diagram of thrust

# type_structure = 'dome'
# type_formdiagram = 'radial_fd'
# discretisation = [20, 16]

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# csv_file = os.path.join(folder, title + '_data.csv')
# thicknesses, solutions = open_csv_row(csv_file, cut_last=False)
# thickness_over_radius = [[], []]
# for el in thicknesses[0]:
#     thickness_over_radius[0].append(el/radius)
#     thickness_over_radius[1].append(el/radius)

# img_graph = None
# diagram_of_thrust(thickness_over_radius, solutions, save=img_graph, fill=True, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=[1.5, 2.0, 2.5, 3.0]).show()

thickness_all = {'cross': [], 'fan': []}
sol_all = {'cross': [], 'fan': []}

GSF_ticks = [1.2, 1.4, 1.6, 1.8, 2.0]

xy_limits = [[0.10, 0.05], [100, 60]]

csv_cross = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/cross_fd.csv'
csv_fan = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/fan_fd.csv'

img_cross = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/individual_cross_fd.pdf'
img_fan = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/individual_fan_fd.pdf'

thickness_all['cross'], sol_all['cross'] = import_csv(csv_cross)

img_graph = img_cross
diagram_of_thrust(thickness_all['cross'], sol_all['cross'], save=img_graph, fill=True, xy_limits=xy_limits, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks).show()

thickness_all['fan'], sol_all['fan'] = import_csv(csv_fan)

thickness_all['fan'][0][-1] -= 0.001
thickness_all['fan'][1][-1] -= 0.001

img_graph = img_fan
diagram_of_thrust(thickness_all['fan'], sol_all['fan'], save=img_graph, fill=True, xy_limits=xy_limits, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks).show()

legends = {'cross': 'cross diagram', 'fan': 'fan diagram'}
# legends = {'cross_fd': [r'orthogonal | $\beta=0^{\circ}$', r'orthogonal | $\beta=20^{\circ}$', r'orthogonal | $\beta=40^{\circ}$'], 'fan_fd': [r'fan-like | $\beta=0^{\circ}$', r'fan-like | $\beta=20^{\circ}$', r'fan-like | $\beta=40^{\circ}$']}
# colors = {'cross_fd': ['C0', 'C0', 'C0'], 'fan_fd': ['C1', 'C1', 'C1']}  # These are C0 and C1 in HEX.
# colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.
colors = {'cross': '#1F77B4', 'fan': '#DA6600'}

thks_plot = []
sols_plot = []
legs_plot = []
col_pts = []
for type_formdiagram in ['cross']:
    thks_plot.append(thickness_all[type_formdiagram])
    sols_plot.append(sol_all[type_formdiagram])
    legs_plot.append(legends[type_formdiagram])
    col_pts.append(colors[type_formdiagram])

img_cross = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/ensemble_cross_fd.pdf'
img_all = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/diagram_thrust/ensemble.pdf'

img_graph = img_cross
diagram_of_multiple_thrust(thks_plot, sols_plot, legs_plot, save=img_graph, limit_state=True, fill=True, xy_limits=xy_limits, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks, colors=col_pts).show()

thks_plot = []
sols_plot = []
legs_plot = []
col_pts = []
for type_formdiagram in ['cross', 'fan']:
    thks_plot.append(thickness_all[type_formdiagram])
    sols_plot.append(sol_all[type_formdiagram])
    legs_plot.append(legends[type_formdiagram])
    col_pts.append(colors[type_formdiagram])

img_graph = img_all
diagram_of_multiple_thrust(thks_plot, sols_plot, legs_plot, save=img_graph, limit_state=True, fill=True, xy_limits=xy_limits, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks, colors=col_pts).show()

# # ------------- plot of the dome diagrams

file_csv1 = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/diagram_radial_spaced.csv'
file_csv2 = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/diagram_flower.csv'

thickness_all = {'radial': [], 'spiral': []}
sol_all = {'radial': [], 'spiral': []}

xy_limits = [[0.10, 0.03], [70, 10]]
GSF_ticks = [1.5, 2.0, 2.5, 3.0]

thickness_all['radial'], sol_all['radial'] = import_csv(file_csv1)

img_graph = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/diagram_radial_spaced.pdf'
diagram_of_thrust(thickness_all['radial'], sol_all['radial'], save=img_graph, fill=True, xy_limits=xy_limits, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks).show()


thickness_all['spiral'], sol_all['spiral'] = import_csv(file_csv2)

img_graph = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/diagram_spiral.pdf'
diagram_of_thrust(thickness_all['spiral'], sol_all['spiral'], save=img_graph, fill=True, xy_limits=xy_limits, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=GSF_ticks).show()

img_all = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/ensemble.pdf'

legends = {'radial': 'radial diagram', 'spiral': 'spiral diagram'}
colors = {'radial': '#42DF89', 'spiral': '#C70039'}

thks_plot = []
sols_plot = []
legs_plot = []
col_pts = []
for type_formdiagram in ['radial', 'spiral']:
    thks_plot.append(thickness_all[type_formdiagram])
    sols_plot.append(sol_all[type_formdiagram])
    legs_plot.append(legends[type_formdiagram])
    col_pts.append(colors[type_formdiagram])

img_graph = img_all
diagram_of_multiple_thrust(thks_plot, sols_plot, legs_plot, save=img_graph, limit_state=True, fill=True, x_label=r'thickness-over-radius ($t/R$)', xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=col_pts).show()

# If wanted to superimpose with the CAS result

# better_csv = '/Users/mricardo/compas_dev/me/min_thk/dome/radial_fd/min_max/dome_radial_fd_discr_[20, 16]_data.csv'

# from compas_tno.plotters import open_csv_row

# thickness_all['radial'], sol_all['radial'] = open_csv_row(better_csv, cut_last=False)
# thickness_over_radius = [[],[]]
# for el in thickness_all['radial'][0]:
#     thickness_over_radius[0].append(el/radius)
#     thickness_over_radius[1].append(el/radius)

# thickness_all['radial'] = thickness_over_radius

# img_graph =  None
# diagram_of_thrust(thickness_all['radial'], sol_all['radial'], save=img_graph, fill=True, limit_state=False, x_label=r'thickness-over-radius ($t/R$)', GSF_ticks=[1.5, 2.0, 2.5, 3.0]).show()

# thks_plot = []
# sols_plot = []
# legs_plot = []
# col_pts = []
# for type_formdiagram in ['radial', 'spiral']:
#     thks_plot.append(thickness_all[type_formdiagram])
#     sols_plot.append(sol_all[type_formdiagram])
#     legs_plot.append(legends[type_formdiagram])
#     col_pts.append(colors[type_formdiagram])

# img_graph = '/Users/mricardo/compas_dev/me/minmax/dome/diagram_thrust/ensemble_diagram_not_spaced.pdf'
# diagram_of_multiple_thrust(thks_plot, sols_plot, legs_plot, save=img_graph, limit_state=True, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=col_pts).show()
