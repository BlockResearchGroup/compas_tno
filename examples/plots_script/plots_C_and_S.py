import os
import matplotlib.pyplot as plt
import numpy as np
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import diagram_of_thrust

# Plot for graph of curves and discretisations
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


type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [20, 16]

folder = os.path.join('/Users/mricardo/compas_dev/me', 'min_thk', type_structure, type_formdiagram, 'min_max')
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

csv_file = os.path.join(folder, title + '_data.csv')
thicknesses, solutions = open_csv_row(csv_file, cut_last=False)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False, GSF_ticks=[1.5, 2.0, 2.5, 3.0]).show()
