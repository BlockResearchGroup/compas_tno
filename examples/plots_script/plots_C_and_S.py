import os
import matplotlib.pyplot as plt
import numpy as np
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import diagram_of_thrust
import csv

size_axis_label = 12
size_axis_data = 10
size_legend = 10

# # Plot for graph of curves and discretisations
# series = []
# x = []
# y = []
# styles = []
# colours = []

# radius = 5.0
# bench = 0.042
# discr = [4, 8, 12, 16, 20, 24]

# # meridians = 12
# series.append(r'$n_\mathrm{M} = 12$')
# y.append(np.array([0.06670, 0.17474, 0.19103, 0.19473, 0.20282, 0.20311])/radius)

# # meridians = 16
# series.append(r'$n_\mathrm{M} = 16$')
# y.append(np.array([0.06922, 0.17687, 0.19266, 0.19669, 0.20455, 0.20464])/radius)

# # meridians = 20
# series.append(r'$n_\mathrm{M} = 20$')
# y.append(np.array([0.07037, 0.17786, 0.19341, 0.19759, 0.20534, 0.20534])/radius)

# # meridians = 24
# series.append(r'$n_\mathrm{M} = 24$')
# y.append(np.array([0.07100, 0.17839, 0.19382, 0.19808, 0.20577, 0.20572])/radius)

# # theoretical
# series.append('benchmark')
# bench_x = [4, 24]
# y.append([bench, bench])

# # for i in range(len(series)):
# #     plt.plot(x[i], y[i], 'o-', label=series[i])

# size_axis_label = 12
# size_axis_data = 10
# size_legend = 10

# fig = plt.figure(figsize=(10, 5))  # try 12, 4
# plt.plot(discr, y[0], 'o-', label=series[0])
# plt.plot(discr, y[1], 'o-', label=series[1])
# plt.plot(discr, y[2], 'o-', label=series[2])
# plt.plot(discr, y[3], 'o-', label=series[3])
# plt.plot(bench_x, y[4], color='black', linestyle='dashed', label=series[4])
# plt.legend(fontsize=size_legend)
# ax = plt.axes()
# ax.set_xlabel(r'number of parallels $(n_\mathrm{P})$', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness over radius $(t/R)$', size=size_axis_label, labelpad=8)
# ax.annotate(r'benchmark $(t/R) = 0.042$', (sum(bench_x)/2, 0.042 + 0.001), textcoords="offset points", xytext=(sum(bench_x)/2, 0.042), ha='center')  # size =
# ax.set_xlim(4-1, 24+1)
# ax.set_ylim(0, 0.05)
# ax.set_xticks([4, 8, 12, 16, 20, 24])
# plt.show()

# ------------------- Plot of Dome

# type_structure = 'dome'
# type_formdiagram = 'radial_fd'
# discretisation = [20, 16]

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

# csv_file = os.path.join(folder, title + '_data.csv')
# thicknesses, solutions = open_csv_row(csv_file, cut_last=False)

# img_graph = os.path.join(folder, title + '_diagram.pdf')
# diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False, GSF_ticks=[1.5, 2.0, 2.5, 3.0]).show()

# ------------------- Plot of Amiens solution

# type_structure = 'crossvault'
# type_formdiagram = 'fan_fd'
# discretisation = 14
# file_name = 'amiens_internet'

# folder = os.path.join('/Users/mricardo/compas_dev/me', 'max_n', file_name, type_structure, type_formdiagram, 'min_max')
# os.makedirs(folder, exist_ok=True)
# title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_offset-method'

# csv_file = os.path.join(folder, title + '_data.csv')
# thicknesses, solutions = open_csv_row(csv_file, cut_last=False)

# img_graph = os.path.join(folder, title + '_diagram.pdf')
# diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False, GSF_ticks=[1.1, 1.2, 1.3, 1.4, 1.5, 1.6]).show()

# ------------------- Plot of Crossvault graph

# path = '/Users/mricardo/compas_dev/me/min_thk/crossvault/study_crossvault_minthk.csv'

# discr = []
# cross_fd = [[], [], []]
# fan_fd = [[], [], []]
# radius = 5.0
# series = ['Analytical min thickness', 'Offset min thickness', 'Correction on normals']
# series = ['Minimum thickness - cross_fd', 'Minimum thickness - fan_fd']

# with open(path) as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line = 0
#     for row in csv_reader:
#         if line > 0:
#             discr.append(int(row[0]))
#             for i in range(3):
#                 cross_fd[i].append(float(row[1 + 2*i])/radius)
#                 fan_fd[i].append(float(row[2 + 2*i])/radius)
#         line += 1

# print(cross_fd)
# print(fan_fd)
# lines = []
# fig = plt.figure(figsize=(12, 5))  # try 12, 4
# ax = plt.subplot(111)
# lines += ax.plot(discr, cross_fd[0], 'o-', label=series[0])#, color='royalblue')
# # lines += ax.plot(discr, cross_fd[1], 's--', label=series[1])#, color='royalblue')
# # lines += ax.plot(discr, cross_fd[2], 'x--', label=series[2])#, color='royalblue')

# lines += ax.plot(discr, fan_fd[0], 'o-', label=series[1])#, color='mediumpurple')
# # lines += ax.plot(discr, fan_fd[1], 's--', label=series[1])#, color='mediumpurple')
# # lines += ax.plot(discr, fan_fd[2], 'x--', label=series[2])#, color='mediumpurple')

# # plt.plot(bench_x, y[4], color='black', linestyle='dashed', label=series[4])
# # ax.legend(lines[:2], ['line A', 'line B'],
# #           loc='upper right', frameon=False)


# from matplotlib.legend import Legend

# # ax = plt.axes()

# ax.legend()

# # ax.legend(lines[:3], series, loc='upper right', fontsize=size_legend, title='Orthogonal Diagram:', bbox_to_anchor=(1.35, 0.8))

# # leg = Legend(ax, lines[3:], series, frameon=False, loc='lower right', fontsize=size_legend, title='Fan-like Diagram:', bbox_to_anchor=(1.35, 0.25))
# # ax.add_artist(leg)

# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

# # leg = Legend(ax, lines[2:], ['line C', 'line D'],
# #              loc='lower right', frameon=False)

# ax.set_xlabel(r'discretisation $(n)$', size=size_axis_label, labelpad=8)
# ax.set_ylabel(r'thickness over radius $(t/R)$', size=size_axis_label, labelpad=8)
# # ax.annotate(r'benchmark $(t/R) = 0.042$', (sum(bench_x)/2, 0.042 + 0.001), textcoords="offset points", xytext=(sum(bench_x)/2, 0.042), ha='center')  # size =
# ax.set_xlim(discr[0]-1, discr[-1]+1)
# ax.set_ylim(0, 0.10)
# ax.set_xticks(discr)
# plt.show()

# --------------------------

path = '/Users/mricardo/compas_dev/me/min_thk/crossvault/study_crossvault_A_D=14_minthk.csv'

A = []
cross_fd = [[], [], [], []]
fan_fd = [[], [], [], []]
radius = 5.0
series = ['analytical', 'from middle', 'from UB/LB', 'from LB']

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 0:
            A.append(float(row[0]))
            for i in range(4):
                cross_fd[i].append(float(row[1 + 2*i])/radius)
                fan_fd[i].append(float(row[2 + 2*i])/radius)
        line += 1

lines = []
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)
lines += ax.plot(A, cross_fd[0], 'o-', label=series[0], color='C0')
lines += ax.plot(A, cross_fd[1], 's--', label=series[1], color='C2')
lines += ax.plot(A, cross_fd[2], 'x--', label=series[2], color='C3')
lines += ax.plot(A, cross_fd[3], '.--', label=series[3], color='C4')

ax.legend()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

ax.set_xlabel(r'reduction coefficient $(A)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'thickness over radius $(t/R)$', size=size_axis_label, labelpad=8)
ax.set_xlim(A[0]-0.025, A[-1]+0.025)
ax.set_ylim(0, 0.18)
ax.set_xticks(A)
plt.show()


lines = []
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)
lines += ax.plot(A, fan_fd[0], 'o-', label=series[0], color='C1')
lines += ax.plot(A, fan_fd[1], 's--', label=series[1], color='C2')
lines += ax.plot(A, fan_fd[2], 'x--', label=series[2], color='C3')
lines += ax.plot(A, fan_fd[3], '.--', label=series[3], color='C4')

ax.legend()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

ax.set_xlabel(r'reduction coefficient $(A)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'thickness over radius $(t/R)$', size=size_axis_label, labelpad=8)
ax.set_xlim(A[0]-0.025, A[-1]+0.025)
ax.set_ylim(0, 0.18)
ax.set_xticks(A)
plt.show()
