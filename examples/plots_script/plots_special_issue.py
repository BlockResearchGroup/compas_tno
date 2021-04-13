import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_multiple_thrust
import csv

size_axis_label = 12
size_axis_data = 10
size_legend = 10

# ------------------- Plot of FAN_FD graph comparison

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/pointed_crossvault_fan_fd_spr_angles.csv'

fan_min_results = [[], [], [], [], []]
fan_min_x = [[], [], [], [], []]

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 1:
            for i in range(5):
                fan_min_results[i].append(float(row[1 + 2*i]))
                fan_min_x[i].append(float(row[0 + 2*i]))
        line += 1


labels =[r'$\beta=0^{\circ}$', r'$\beta=10^{\circ}$', r'$\beta=20^{\circ}$', r'$\beta=30^{\circ}$', r'$\beta=40^{\circ}$']
ticks = ['o-', 's--', 'x--', ]
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)
for i in range(5):
    ax.plot(fan_min_x[i], fan_min_results[i], label=labels[i])

ax.set_xlabel(r'radius over length $(r/l_0)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'minimum thickness over span $(t_{min}/s)$', size=size_axis_label, labelpad=8)
ax.set_ylim(0, 0.06)
# ax.set_xticks([0.50, 0.75, 1.0, 1.25, 1.5])
ax.set_xticks([0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.10, 1.20, 1.30, 1.40, 1.50])
ax.legend()
plt.show()

# print('\nFAN 3D')
# print(fan_min_results)
# print(fan_min_x)

# degrees = []
# i = 0
# for serie in fan_min_x:
#     n = len(serie)
#     degrees.append([i]*n)
#     i += 10

# Xfan = np.array(fan_min_x)
# Yfan = np.array(degrees)
# Zfan = np.array(fan_min_results)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf1 = ax.plot_surface(Xfan, Yfan, Zfan, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf1, shrink=0.5, aspect=5)
# plt.show()

# # ------------------- Plot of CROSS_FD graph comparison

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/pointed_crossvault_cross_fd_spr_angles.csv'

cross_min_results = [[], [], [], [], []]
cross_min_x = [[], [], [], [], []]

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 1:
            for i in range(5):
                if row[1 + 2*i] == '':
                    continue
                cross_min_results[i].append(float(row[1 + 2*i]))
                cross_min_x[i].append(float(row[0 + 2*i]))
        line += 1

labels =[r'$\beta=0^{\circ}$', r'$\beta=10^{\circ}$', r'$\beta=20^{\circ}$', r'$\beta=30^{\circ}$', r'$\beta=40^{\circ}$']
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)
for i in range(5):
    ax.plot(cross_min_x[i], cross_min_results[i], label=labels[i])

ax.set_xlabel(r'radius over length $(r/l_0)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'minimum thickness over span $(t_{min}/s)$', size=size_axis_label, labelpad=8)
ax.set_ylim(0, 0.06)
ax.set_xticks([0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.10, 1.20, 1.30, 1.40, 1.50])
ax.legend()
plt.show()


# print('\CROSS 3D')
# print(cross_min_results)
# print(cross_min_x)

# for j in range(3):
#     cross_min_results[4].append(0.08)
#     cross_min_x[4].append(cross_min_x[3][j + 19])

# degrees = []
# i = 0
# for serie in cross_min_x:
#     n = len(serie)
#     degrees.append([i]*n)
#     i += 10

# Xcross = np.array(cross_min_x)
# Ycross = np.array(degrees)
# Zcross = np.array(cross_min_results)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf2 = ax.plot_surface(Xcross, Ycross, Zcross, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf2, shrink=0.5, aspect=5)
# plt.show()

# # # ------------------- Plot Fan and Cross 20 DEG

labels =[r'fan-like | $\beta=20^{\circ}$', r'orthogonal | $\beta=20^{\circ}$']
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)

ax.plot(fan_min_x[2], fan_min_results[2], label=labels[0])
ax.plot(cross_min_x[2], cross_min_results[2], label=labels[1])

ax.set_xlabel(r'radius over length $(r/l_0)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'minimum thickness over span $(t_{min}/s)$', size=size_axis_label, labelpad=8)
ax.set_ylim(0, 0.06)
ax.set_xticks([0.50, 0.75, 1.0, 1.25, 1.5])
ax.legend()
plt.show()

# # # ------------------- Plot Fan + Cross + Cross Braced for 20 DEG

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/pointed_crossvault_crossbraced_20deg.csv'

labels.append(r'curved | $\beta=20^{\circ}$')
fig = plt.figure(figsize=(12, 5))  # try 12, 4
ax = plt.subplot(111)

braced_x = [[]]
braced_results = [[]]

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 1:
            for i in range(1):
                if row[1 + 2*i] == '':
                    continue
                braced_results[i].append(float(row[1 + 2*i]))
                braced_x[i].append(float(row[0 + 2*i]))
        line += 1

ax.plot(fan_min_x[2], fan_min_results[2], label=labels[0])
ax.plot(cross_min_x[2], cross_min_results[2], label=labels[1])
ax.plot(braced_x[0], braced_results[0], label=labels[2])

ax.set_xlabel(r'radius over length $(r/l_0)$', size=size_axis_label, labelpad=8)
ax.set_ylabel(r'minimum thickness over span $(t_{min}/s)$', size=size_axis_label, labelpad=8)
ax.set_ylim(0, 0.06)
ax.set_xticks([0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.10, 1.20, 1.30, 1.40, 1.50])
ax.legend()
plt.savefig('/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/MDPI/beta20diagrams.pdf')
plt.show()


# ------------- plot of the diagrams of thrust

# xy_limits = [[0.050, 0.001], [100, 50]]
# GSF_ticks = [2.0, 3.0, 4.0, 5.0]
# deg = 20
# Rs = [6.1147, 7.08, 7.9422]
# # legends = {'cross_fd': [r'orthogonal | $\beta=20$'], 'fan_fd': [r'fan-like | $\beta=20$'], 'topology-crossbraced': [r'braced | $\beta=20$']}
# # colors = {'cross_fd': ['C0'], 'fan_fd': ['C1'], 'topology-crossbraced': ['C2']}  # These are C0 and C1 in HEX.
# colors = ['C0', 'C1', 'C2']
# legends = [r'fan-like | $\beta=20$', r'orthogonal | $\beta=20$', r'curved | $\beta=20$']
# # colors = {'cross_fd': ['#419EDE', '#1F77B4', '#144C73'], 'fan_fd': ['#FFA85B', '#FF7F0E', '#C15A00']}  # These are C0 and C1 in HEX.
# # colors = {'cross_fd': ['#1FB4A7', '#1F77B4', '#1F2DB4'], 'fan_fd': ['#FFA85B', '#DA6600', '#FF0E16']}  # These are C0 and C1 in HEX.

# type_structure = 'pointed_crossvault'
# discretisation = 14
# option = ['A', 'B', 'C']

# i_ = 0
# for R in Rs:
#     thicknesses_all = []
#     solutions_all = []
#     for type_formdiagram in ['fan_fd', 'cross_fd', 'topology-crossbraced']:
#         folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram, 'R='+str(R), 'min_thk', 'deg='+str(deg))
#         # os.makedirs(folder, exist_ok=True)
#         title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
#         forms_address = os.path.join(folder, title)
#         csv_file = os.path.join(folder, title + '_data.csv')
#         thicknesses, solutions = open_csv_row(csv_file, cut_last=True)
#         img_graph = None
#         adim_thk = [[], []]
#         print(thicknesses)
#         for i in range(len(thicknesses)):
#             for j in range(len(thicknesses[i])):
#                 adim_thk[i].append(thicknesses[i][j]/10.0)
#         print(adim_thk)
#         # diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, limit_state=False).show()
#         thicknesses_all.append(adim_thk)
#         solutions_all.append(solutions)
#         print(type_formdiagram, deg, -solutions[1][0]/solutions[0][0])

#     folder_main = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, 'MDPI')
#     title_main = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
#     img_graph = os.path.join(folder_main, title_main + '_diagram_MDPI_DIAGRAM_'+option[i_]+'.pdf')
#     diagram_of_multiple_thrust(thicknesses_all, solutions_all, legends, save=img_graph, fill=True, xy_limits=xy_limits, GSF_ticks=GSF_ticks, colors=colors).show()
#     i_ += 1
