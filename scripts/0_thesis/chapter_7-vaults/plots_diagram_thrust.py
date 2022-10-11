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
        markers = ['', '', '', '', '', '', '', '', '', '', '4', '8']

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

    return fig, ax

size_axis_label = 12
size_axis_data = 10
size_legend = 10

xy_limits = [[0.050, 0.001], [100, 50]]
GSF_ticks = [2.0, 3.0, 4.0, 5.0]
deg = 20
Rs = [6.1147, 7.08, 7.9422]
# legends = {'cross_fd': [r'orthogonal | $\beta=20$'], 'fan_fd': [r'fan-like | $\beta=20$'], 'topology-crossbraced': [r'braced | $\beta=20$']}
# colors = {'cross_fd': ['C0'], 'fan_fd': ['C1'], 'topology-crossbraced': ['C2']}  # These are C0 and C1 in HEX.
colors = ['C0', 'C1', 'C2']
legends = [r'fan-like | $\beta=20$', r'orthogonal | $\beta=20$', r'curved | $\beta=20$']
# colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

type_structure = 'pointed_crossvault'
discretisation = 14
option = ['A', 'B', 'C']

param_thks = {}
param_mins = {}
param_maxs = {}

curv_thks = {}
curv_mins = {}
curv_maxs = {}

x_add = {}
y_add = {}

# adding parameteric envelope by hand

param_thks['6.1147'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.0138581687892141]
param_mins['6.1147'] = [0.696302181430962, 0.712190170707322, 0.728715423421778, 0.746097517580299, 0.765045434280777, 0.785257690737793, 0.806267884905145, 0.83196507678578, 0.836965223876471]
param_maxs['6.1147'] = [-1.00891403607344, -0.983829961062245, -0.959678524921317, -0.936004852009428, -0.91283679868865, -0.890979621242108, -0.869560153457901, -0.846864654473703, -0.836965223876471]

param_thks['7.08'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.0180566482731426]
param_mins['7.08'] = [0.622129590045542, 0.635886376449303, 0.650113852193332, 0.665025993190104, 0.680726634985757, 0.697321111094259, 0.720838836706226, 0.754437107668423]
param_maxs['7.08'] = [-0.874918142130951, -0.855463033424789, -0.836654226378205, -0.818451721428526, -0.800817948896633, -0.783702346353567, -0.766800192517039, -0.754437107668423]

param_thks['7.9422'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.0210616141712768]
param_mins['7.9422'] = [0.574844815391402, 0.587220045811639, 0.600260450626909, 0.613819594404171, 0.628155133426835, 0.653045545746104, 0.704111404462089]
param_maxs['7.9422'] = [-0.793821830005611, -0.777418630549045, -0.761524986455031, -0.746111584081397, -0.731150795765791, -0.716616504228468, -0.704111404462089]

legends.append(r'parametric env. | $\beta=20$')
colors.append('C3')

# adding envelope of curved diagrams by hand

curv_thks['6.1147'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.0138581661416226]
curv_mins['6.1147'] = [0.696302181430962, 0.712190170707322, 0.728715423421778, 0.746097517580299, 0.765045434280777, 0.785618227177949, 0.807334674761319, 0.83196507678578, 0.839221614598699]
curv_maxs['6.1147'] = [-1.08311049024263, -1.05618171549355, -1.03025416535101, -1.00525899957417, -0.981131952900054, -0.957812961802496, -0.93159719563686, -0.859265971676013, -0.839221614598699]

curv_thks['7.08'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.0203592089760973]
curv_mins['7.08'] = [0.623080862340747, 0.637756507956512, 0.65303731323751, 0.668965085538148, 0.685585354937527, 0.702947815919594, 0.780113196617196]
curv_maxs['7.08'] = [-0.939456990462269, -0.918566775789761, -0.898370540407471, -0.878825351599485, -0.859890778089617, -0.841528852891464, -0.780113196617196]

curv_thks['7.9422'] = [0.05, 0.045, 0.04, 0.035, 0.03, 0.0270472443995507]
curv_mins['7.9422'] = [0.578634829095247, 0.591704212836187, 0.605291881945596, 0.619432139493077, 0.640535925881035, 0.696513170372252]
curv_maxs['7.9422'] = [-0.8524907198411, -0.83487520924929, -0.817806904934335, -0.801254328788241, -0.754650279026647, -0.696513170372252]


x_add['6.1147'] = [0.0139, 0.0142, 0.0146, 0.0151, 0.0155, 0.0159, 0.0164, 0.0169]
y_add['6.1147'] = [0.8392, 0.8475, 0.8557, 0.8641, 0.8725, 0.8814, 0.8910, 0.9008]

x_add['7.08'] = []
y_add['7.08'] = []

x_add['7.9422'] = []
y_add['7.9422'] = []

legends.append(r'curved env. | $\beta=20$')
colors.append('C4')

i_ = 0
for R in Rs:
    thicknesses_all = []
    solutions_all = []
    for type_formdiagram in ['fan_fd', 'cross_fd', 'topology-crossbraced']:
        folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(R), 'min_thk', 'deg='+str(deg))
        # os.makedirs(folder, exist_ok=True)
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
        forms_address = os.path.join(folder, title)
        csv_file = os.path.join(folder, title + '_data.csv')
        print('File with data:', csv_file)
        thicknesses, solutions = open_csv_row(csv_file, cut_last=True)
        img_graph = None
        adim_thk = [[], []]
        for i in range(len(thicknesses)):
            for j in range(len(thicknesses[i])):
                adim_thk[i].append(thicknesses[i][j]/10.0)
        print(adim_thk)
        # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
        thicknesses_all.append(adim_thk)
        solutions_all.append(solutions)
        print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])

    thicknesses_all.append([param_thks[str(R)], param_thks[str(R)]])
    solutions_all.append([param_mins[str(R)], param_maxs[str(R)]])

    thicknesses_all.append([curv_thks[str(R)], curv_thks[str(R)]])
    solutions_all.append([curv_mins[str(R)], curv_maxs[str(R)]])

    folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, 'MDPI')
    title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
    img_graph = os.path.join(folder_main, title_main + 'NOSAVE'+option[i_]+'.pdf')

    print('thicknesses', thicknesses_all)
    print('solutionss', solutions_all)

    fig, ax = diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors)

    ax.plot(x_add[str(R)], y_add[str(R)])
    plt.show()

    i_ += 1
