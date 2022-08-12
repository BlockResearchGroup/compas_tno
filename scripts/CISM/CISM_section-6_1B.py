from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.algorithms import compute_reactions
from compas_tno.problems import initialise_form
from compas_tno.problems import adapt_problem_to_fixed_diagram
from compas_tno.utilities import apply_envelope_from_shape
from compas.colors import Color

from matplotlib.ticker import FormatStrFormatter
from numpy import arange
from numpy import array
from numpy import append
from compas.geometry import intersection_line_line_xy
import matplotlib.pyplot as plt

def diagram_of_thrust(thicknesses, solutions, limit_state=True, fill=True, xy_limits=None, x_label=None, show_legend=True, show_GSFscale=False, GSF_ticks=False, save=False):
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

    interval_x = 0.005
    interval_y = 20
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


delta = 1.0
span = 10.0
xspan = yspan = [0.0, span]
xspan_vault = yspan_vault = [- delta, span + delta]

# from compas_tno.utilities import diagram_of_thrust
import compas_tno

discr = 16
path = compas_tno.get('')

thicknesses = []
solutions = []

i = 0
for objective in ['min', 'max']:
    thicknesses.append([])
    solutions.append([])
    for thk in [0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.067]:

        if thk == 0.067:
            form_path = '//Users/mricardo/compas_dev/compas_tno/data/CISM/form-t-16.json'
        else:
            form_path = path + '/CISM/form-' + objective + '-' + str(discr) + '-thk-' + str(thk) + '.json'
        print(form_path)
        form = FormDiagram.from_json(form_path)
        thrust = form.thrust()
        swt = form.lumped_swt()
        T_over_W = abs(thrust/swt)

        thicknesses[i].append(thk/10)
        solutions[i].append(T_over_W)

    i += 1

print(thicknesses)
print(solutions)

plotter = diagram_of_thrust(thicknesses, solutions)
plotter.show()

