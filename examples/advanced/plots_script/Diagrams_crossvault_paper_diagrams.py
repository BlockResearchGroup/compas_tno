import os
from compas_tno.plotters import open_csv
from compas_tno.plotters import interpolate_min_thk
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import prune_data
from compas_tno.plotters import diagram_of_multiple_thrust
import matplotlib.pyplot as plt
from numpy import array

# ['o', '^', 's', 'D', 'x', '1', '2', '3', 'v', 'p', '4', '8']


def GSF_plot(Parameters, GSF, legends, markers, colors, save=None):

    fig = plt.figure(figsize=[8, 6])
    ax = fig.add_subplot()

    for i in range(len(Parameters)):
        print(Parameters[i])
        print(GSF[i])
        ax.plot(Parameters[i], GSF[i], markers[i], ls='-', color=colors[i], label=legends[i])

    ax.set_xlabel('Radius', weight='bold')
    ax.set_ylabel('GSF', weight='bold')
    ax.legend()
    if save:
        plt.savefig(save)
    plt.show()

    return


xs = []
mins = []
maxs = []
legends = []
colors = []

thk = 0.6
span = 10.0
type_structures = ['pointed_crossvault']
type_mesh = 'mix'
type_formdiagrams = ['cross_fd', 'fan_fd', 'topology-mix']  # ['cross_fd', 'fan_fd', 'cross_fd', 'fan_fd', 'topology-mix']
sags = [False, False, False] # [False, False, 50.0, 10.0, False]
# type_formdiagrams = ['topology-mix']
hcs = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]  # [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]
Rs = [5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
colors_list = ['#D2B4DE', '#b4c0de', '#a1eb85']  # ['#D2B4DE', '#b4c0de', 'orange', 'blue', '#a1eb85']
# colors_list = ['#a1eb85']
discretisation = 10
markers = ['o', '^', 's', 'D', 'x', '1', '2', '3']*10
xy_limits = [[0.60, 0.20], [110, 30]]
legend = False

GSF = []
RS = []
GSF_Legends = []
i = 0

for type_structure in type_structures:
    for type_formdiagram in type_formdiagrams:
        R = 5
        gsf_pattern = []
        rs_pattern = []
        for hc in hcs:
            if type_formdiagram == 'topology-mix':
                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
                title = type_structure + '_' + type_formdiagram + '_Mesh-' + type_mesh + '_discr_' + str(100) + 'smooth_'
                csv_file = os.path.join(folder, title + '_data.csv')
            else:
                folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
                title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
                if sags[i]:
                    title = title + 'sag_' + str(sags[i])
                csv_file = os.path.join(folder, title + '_data.csv')
            print(csv_file)
            sizes, solutions = open_csv_row(csv_file, printout=False)
            sizes, solutions = prune_data(sizes, solutions)
            size_parameters = sizes[0]
            solutions_min = solutions[0]
            solutions_max = solutions[1]
            xs.append(size_parameters)
            mins.append(solutions_min)
            maxs.append(solutions_max)
            min_thk = interpolate_min_thk(sizes, solutions)
            gsf_pattern.append(thk/min_thk)
            rs_pattern.append(R)
            colors.append(colors_list[i])
            legends.append('Crossvault R='+str(R)+'_'+type_formdiagram)
            R += 0.5
        GSF.append(gsf_pattern)
        RS.append(rs_pattern)
        if sags[i]:
            GSF_Legends.append(type_formdiagram + '_sag_' + str(sags[i]))
        else:
            GSF_Legends.append(type_formdiagram)
        i += 1

print('Check on the parameters size')
print(len(xs))
print(len(mins))
print(len(maxs))
print(len(colors))
print(len(legends))

print(GSF)
print(RS)
GSF_file = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, 'GSF_plot' + str(type_formdiagrams) + '.pdf')
GSF_plot(RS, GSF, GSF_Legends, markers, colors_list, save=GSF_file)

diagram_file = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, 'diagram_GSF_'+ str(type_formdiagrams) +'.pdf')
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=legend, xy_limits=xy_limits, markers=markers).show()
