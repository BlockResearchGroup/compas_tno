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


def diagram_of_thrust(thicknesses, solutions, limit_state=True, fill=True, xy_limits=None, x_label=None, show_legend=True, show_GSFscale=False, GSF_ticks=False, save=False, interval_x = 0.05, interval_y = 20):
    """Plot a diagram of Thrusts based on the collected data from (n) points.

    Parameters
    ----------
    thicknesses : list of lists [[n],[m]]
        Points with discretised solutions for minimum/maximum thrust .
    solutions : list of lists [[n],[m]]
        Adimensional thrust over weight for minimum/maximum thrust.
    limit_state : bool
        If yes, it interpolates the limit_state based on the two last solutions.
    fill : bool
        If yes, it fills bewteen the curves of min and max.
    xy_limits : list of lists
        Freezes the axis of the graphs.
    x_label : str, optional
        Whether or not labels should be activated in the x-axis, by default None
    show_legend : bool, optional
        Whether or not legend should be shown, by default True
    show_GSFscale : bool, optional
        Whether or not show the GSF scale, by default False
    GSF_ticks : bool, optional
        List with the ticks for the GSF plot, by default False
    save : str, optional
        Whether plot is to besaved, by default False
    interval_x : float, optional
        Default ticks in the x-axis, by default 0.05
    interval_y : int, optional
        Default ticks in the y-axis, by default 20

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
    fmax = 100.0 * array(max_sol)
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

    fig = plt.figure(figsize=[10, 4])

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
        ax.annotate(str(round(max_x/x_, 1)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    if fill:
        ax.fill(append(xmin_, xmax_[::-1]), append(fmin_, fmax_[::-1]), color="grey", alpha=0.2)

    # ax1 = plt.axes()
    if show_GSFscale:
        ax2 = ax.twiny()
        ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
        ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
        ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)

    if x_label:
        ax.set_xlabel(x_label, size=size_axis_label, labelpad=8)
    else:
        ax.set_xlabel(r'$\mathbf{t/s}$', size=size_axis_label, labelpad=8)
    ax.set_ylabel(r'$\mathbf{T_i/W_i}$ [%]', size=size_axis_label, labelpad=8)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # Check if this is necessary
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    ax.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax.tick_params(axis='both', which='major', labelsize=size_axis_data)

    box = ax.get_position()
    ax.set_position([box.x0*0.7, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    if show_legend:
        ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt


def diagram_of_energy(thicknesses, solutions, limit_state=True, fill=True, xy_limits=None, x_label=None, show_legend=True, show_GSFscale=False, GSF_ticks=False, save=False):
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
    fmax = 100.0 * array(max_sol)
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

    interval_x = 0.05
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

    fig = plt.figure(figsize=[10, 4])

    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xmax, fmax, 'o', ls='-', markersize=5, color='gray', label='inwards energy')
    ax.plot(xmin, fmin, 'o', ls='-', markersize=5, color='orange', label='outwards energy')

    if limit_state:
        extrapolation_max = [fmax[m-1], y_]
        extrapolation_min = [fmin[n-1], y_]
        extrapolation_x_max = [xmax[m-1], x_]
        extrapolation_x_min = [xmin[n-1], x_]
        ax.plot(x_, y_, 'o', ls=' ', markersize=7, color='black', label='limit state')
        ax.plot(extrapolation_x_max, extrapolation_max, color='gray')
        ax.plot(extrapolation_x_min, extrapolation_min, color='orange')
        # ax.annotate(str(round(max_x/x_, 2)), (x_, y_), textcoords="offset points", xytext=(0, 10), ha='center')

    if fill:
        ax.fill(append(xmin_, xmax_[::-1]), append(fmin_, fmax_[::-1]), color="grey", alpha=0.2)

    # ax1 = plt.axes()
    if show_GSFscale:
        ax2 = ax.twiny()
        ax2.set_xticks([0] + [100*(max_x-tck_x)/(max_x-min_x) for tck_x in ticks_x] + [100])
        ax2.set_xticklabels(['1.0'] + [str(round(tck_GSF, 2)) for tck_GSF in ticks_GSF] + [''], size=size_axis_data)
        ax2.set_xlabel('GSF', size=size_axis_label, labelpad=8)

    if x_label:
        ax.set_xlabel(x_label, size=size_axis_label, labelpad=8)
    else:
        ax.set_xlabel(r'$\mathbf{t}$', size=size_axis_label, labelpad=8)
    ax.set_ylabel(r'$\mathbf{E_i/W_i}$ [%]', size=size_axis_label, labelpad=8)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # Check if this is necessary
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlim(max_x, min_x)
    ax.set_ylim(min_y, max_y)
    ax.set_xticks(arange(max_x, min_x - interval_x, -interval_x))
    ax.set_yticks(arange(min_y, max_y + interval_y, interval_y))
    ax.tick_params(axis='both', which='major', labelsize=size_axis_data)

    box = ax.get_position()
    ax.set_position([box.x0*0.7, box.y0*1.5, box.width * 0.90, box.height*0.90])
    ax.grid(color='silver', linestyle='-', linewidth=0.5)
    if show_legend:
        ax.legend(fontsize=size_legend)  # , bbox_to_anchor(0.1, 0.1), ncol=1)

    if save:
        plt.savefig(save)

    return plt



def box_limits(box):
    return [box.x0, box.y0, box.width*0.8, box.height]


size_plots = (8, 5)

size_axis_label = 12
size_axis_data = 12
size_legend = 12

# ## -----------------------
# ## - Data of the anagni vault
# ## -----------------------

# thk = [0.25, 0.23, 0.21, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.07, 0.056]
# min_thrust = [0.753, 0.77, 0.79, 0.81, 0.83, 0.86, 0.88, 0.91, 0.94, 0.98, 1.04]
# max_thrust = [1.248, 1.22, 1.20, 1.18, 1.16, 1.14, 1.12, 1.10, 1.08, 1.06, 1.04]

# xy_limits = [[0.25, 0.0], [130, 70]]
# GSF_ticks = [2.0, 3.0, 4.0, 5.0]

# diagram_of_thrust([thk, thk], [min_thrust, max_thrust], xy_limits=xy_limits, GSF_ticks=GSF_ticks).show()


# Data from the test on the vault under 1-corner displacement

## -----------------------
## Data from the test on the vault under 1-corner displacement
## -----------------------

thk = [0.5000, 0.4500, 0.4000, 0.3500, 0.3355]
thk2 = [thk, thk]

# max_thrust = [1.020, 0.970, 0.923, 0.865, 0.844]
# min_thrust = [0.760, 0.783, 0.808, 0.834, 0.844]

# xy_limits = [[0.50, 0.30], [120, 60]]

# diagram_of_thrust([thk, thk], [min_thrust, max_thrust], xy_limits=xy_limits).show()

# diagram_of_energy([thk, thk], [min_thrust, max_thrust], xy_limits=xy_limits).show()


## -----------------------
## Data from the test on the vault under 1-corner displacement
## -----------------------

xy_limits = [[0.50, 0.30], [50, 0]]

E_Wmin = {}
T_Wmin = {}
E_Wmax = {}
T_Wmax = {}

T_Wmin[0] = [0.852617799596812, 0.854101861253407, 0.840679448707132, 0.835296029760983, 0.842248520137569]
E_Wmin[0] = [0.163883470777335, 0.170225486835804, 0.178333561146456, 0.19853514140168, 0.209532236491124]
T_Wmax[0] = [0.965954411053202, 0.927286869598406, 0.886710858802386, 0.85385423864999, 0.845019374201551]
E_Wmax[0] = [0.289688467140803, 0.270220747029148, 0.25085041463346, 0.221236190719395, 0.212213131041641]

T_Wmin[22.5] = [0.85408081858013, 0.854101861442466, 0.840675934432547, 0.835295562036828, 0.842247344285949]
E_Wmin[22.5] = [0.240540191912474, 0.247035263122189, 0.255682270797365, 0.277569497307334, 0.289122569047829]
T_Wmax[22.5] = [0.966122127756217, 0.919914599217358, 0.886824195661815, 0.850249784602811, 0.845016898636713]
E_Wmax[22.5] = [0.373178831981132, 0.353705601551965, 0.334020159506617, 0.302082886388867, 0.291903089230211]

T_Wmin[45.0] = [0.855284658227807, 0.854101861125143, 0.84068121539555, 0.835328005484981, 0.842244686298617]
E_Wmin[45.0] = [0.280568085425574, 0.286236159775127, 0.294105672605692, 0.314346413421779, 0.324696611369729]
T_Wmax[45.0] = [0.955375851750082, 0.919811130244569, 0.886824787715173, 0.850250567443381, 0.845036051746596]
E_Wmax[45.0] = [0.400244663906881, 0.38357563078471, 0.366338363024731, 0.336960568418875, 0.327153448191514]

T_Wmin[67.5] = [0.93140631585553, 0.90337187629426, 0.875432619738802, 0.835292651688774, 0.842246915688296]
E_Wmin[67.5] = [0.275554421477574, 0.280713307861044, 0.28711710283448, 0.303266937650157, 0.310838537992659]
T_Wmax[67.5] = [0.955373310226573, 0.919909653492312, 0.88682022056708, 0.850251703278628, 0.845009942864741]
E_Wmax[67.5] = [0.366425005733394, 0.355049747352146, 0.342884871492607, 0.320539058354002, 0.312597660172167]

T_Wmin[90.0] = [0.935834854172798, 0.905116214366648, 0.876659378227649, 0.848900168799814, 0.843770854180785]
E_Wmin[90.0] = [0.222740843498805, 0.227111623547614, 0.232393293536222, 0.244681617453372, 0.249542863240051]
T_Wmax[90.0] = [0.931276050843824, 0.905082185162974, 0.876667558019345, 0.850224395200505, 0.843985039405492]
E_Wmax[90.0] = [0.277388399294086, 0.272888376369578, 0.267606706564648, 0.25531838254663, 0.250451712341471]

T_Wmin[-22.5] = [0.818091015126128, 0.851019199474235, 0.840692927018947, 0.835139464961215, 0.842024506585788]
E_Wmin[-22.5] = [0.0621203317399916, 0.0674992455910197, 0.0738351859019236, 0.089275615338085, 0.0980425276226284]
T_Wmax[-22.5] = [0.96395885318127, 0.926631037872626, 0.887382946140621, 0.853171404510166, 0.847109984943211]
E_Wmax[-22.5] = [0.162095655940245, 0.145703098204884, 0.129543228994998, 0.106904742445738, 0.100350842692685]

T_Wmin[-45] = [0.798819298947986, 0.818243382640541, 0.840696621451932, 0.835136921466002, 0.842024766022851]
E_Wmin[-45] = [0.0498515471364692, 0.0461824015017301, 0.041903930063664, 0.03357532095971, 0.0283732762812016]
T_Wmax[-45] = [0.96523767506389, 0.925966636506193, 0.886914978908928, 0.853206842867256, 0.847110022565163]
E_Wmax[-45] = [0.00982525425180188, 0.000996523556246222, 0.0114727945230315, 0.0237019781570626, 0.0267781327297162]

T_Wmin[-67.5] = [0.839585627665333, 0.837073219613384, 0.840655890579557, 0.845904605335619, 0.843149542673727]
E_Wmin[-67.5] = [0.162451084045708, 0.159040683489061, 0.155395260492879, 0.151633733954817, 0.150483414299278]
T_Wmax[-67.5] = [0.985900234587391, 0.933734281546767, 0.890111244146079, 0.85389550324132, 0.847109657841364]
E_Wmax[-67.5] = [0.131981058905596, 0.136279075950859, 0.140913190361772, 0.147624631934958, 0.149783236231981]

for angle in [0, 22.5, 45, 67.5, 90.0]:
    diagram_of_energy([thk, thk], [E_Wmin[angle], E_Wmax[angle]], xy_limits=xy_limits, show_legend=False).show()

print('Next one')

for angle in [-22.5, -45, -67.5]:
    diagram_of_energy([thk, thk], [E_Wmin[angle], E_Wmax[angle]], xy_limits=xy_limits, show_legend=False).show()
