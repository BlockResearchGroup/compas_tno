import os
from compas_tno.plotters import open_csv
from compas_tno.plotters import interpolate_min_thk
from compas_tno.plotters import diagram_of_multiple_thrust
import matplotlib.pyplot as plt
from numpy import array

# ['o', '^', 's', 'D', 'x', '1', '2', '3', 'v', 'p', '4', '8']


def GSF_plot(Parameters, GSF, legends, markers, colors):

    fig = plt.figure(figsize=[8, 6])
    ax = fig.add_subplot()

    for i in range(len(Parameters)):
        print(Parameters[i])
        print(GSF[i])
        ax.plot(Parameters[i], GSF[i], markers[i], ls='-', color=colors[i], label=legends[i])

    ax.set_xlabel('Radius', weight='bold')
    ax.set_ylabel('GSF', weight='bold')
    ax.legend()

    plt.show()

    return


xs = []
mins = []
maxs = []
legends = []
colors = []

thk = 0.5
span = 10.0
type_structures = ['pointed_crossvault']
type_formdiagrams = ['cross_fd']  # , 'fan_fd'
hcs = [5.00]  # , 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]
Rs = [5, 6, 7, 8, 9, 10, 11, 12]
colors_list = ['#D2B4DE', '#b4c0de']
colors_list_shade1 = ['#b27fc7', '#c7a2d6', '#ddc6e6', '#f2e9f5']
colors_list_shade1 = ['#1d0f22', '#4a2758', '#8f4baa', '#b98bcc']
discretisation = 10
markers = ['o', '^', 's', 'D', 'x', '1', '2', '3']*2
markers = ['']*20
xy_limits = [[0.050, 0.02], [110, 60]]
legend = False
hor_combs = [[0.0, 0.0], [1.0, 1.0], [2.5, 2.5], [10.0, 10.0]]

GSF = []
RS = []
GSF_Legends = []

for type_structure in type_structures:
    for type_formdiagram in type_formdiagrams:
        R = 5
        gsf_pattern = []
        rs_pattern = []
        for hc in hcs:
            folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'h='+str(hc))
            i = 0
            for horizontal_force_x, horizontal_force_y in hor_combs:
                if horizontal_force_x == 0.0 and horizontal_force_y == 0.0:
                    title = title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
                    legends.append('No Hor. React')
                else:
                    title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_partial_ratio_' + str(horizontal_force_x) + '_' + str(horizontal_force_y)
                    legends.append('Hor. React='+str([horizontal_force_x, horizontal_force_y])+'_'+type_formdiagram)
                csv_file = os.path.join(folder, title + '_data.csv')
                size_parameters, solutions_min, solutions_max = open_csv(csv_file)
                xs.append(size_parameters)
                mins.append(solutions_min)
                maxs.append(solutions_max)
                min_thk = interpolate_min_thk(size_parameters, solutions_min, solutions_max)
                gsf_pattern.append(thk/min_thk/span)
                rs_pattern.append(R)
                colors.append(colors_list_shade1[i])
                i += 1
            R += 1
        GSF.append(gsf_pattern)
        RS.append(rs_pattern)
        GSF_Legends.append(type_formdiagram)

print('Check on the parameters size')
print(len(xs))
print(len(mins))
print(len(maxs))
print(len(colors))
print(len(legends))

print(GSF)
print(RS)
# GSF_plot(RS, GSF, GSF_Legends, markers, colors_list)

diagram_file = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, 'diagram_hor.pdf')
diagram_of_multiple_thrust(xs, mins, maxs, legends, colors=colors, save=diagram_file, fill=True, show_legend=legend, xy_limits=xy_limits, markers=markers).show()
