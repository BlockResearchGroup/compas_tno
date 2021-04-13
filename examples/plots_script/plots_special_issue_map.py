import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from compas_tno.plotters import open_csv_row
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import diagram_of_multiple_thrust
from matplotlib.ticker import PercentFormatter
import csv

size_axis_label = 12
size_axis_data = 10
size_legend = 10

cmap = 'seismic'

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


print('\nFAN 3D')
print(fan_min_results)
print(fan_min_x)

degrees = []
i = 0
for serie in fan_min_x:
    n = len(serie)
    degrees.append([i]*n)
    i += 10

Xfan = np.array(fan_min_x)
Yfan = np.array(degrees)
Zfan = np.array(fan_min_results)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf1 = ax.plot_surface(Xfan, Yfan, Zfan, cmap=cmap,
#                        linewidth=0, antialiased=False, vmin=0.0, vmax=0.05)
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

print('\CROSS 3D')
print(cross_min_results)
print(cross_min_x)

for j in range(3):
    cross_min_results[4].append(0.08)
    cross_min_x[4].append(cross_min_x[3][j + 19])

degrees = []
i = 0
for serie in cross_min_x:
    n = len(serie)
    degrees.append([i]*n)
    i += 10

Xcross = np.array(cross_min_x)
Ycross = np.array(degrees)
Zcross = np.array(cross_min_results)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf2 = ax.plot_surface(Xcross, Ycross, Zcross, cmap=cmap,
#                        linewidth=0, antialiased=False, vmin=0.0, vmax=0.05)
# fig.colorbar(surf2, shrink=0.5, aspect=5)
# plt.show()

# Load special lines

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/special-lines.csv'

betas = [0, 10, 20, 30, 40]
mincross = []
mincross_x = []
minfan = []
minfan_x = []
minint = []
minint_x = []

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 1:
            i = 0
            mincross.append(float(row[1 + 2*i]))
            mincross_x.append(float(row[0 + 2*i]))
            i = 1
            minfan.append(float(row[1 + 2*i]))
            minfan_x.append(float(row[0 + 2*i]))
            i = 2
            minint.append(float(row[1 + 2*i]))
            minint_x.append(float(row[0 + 2*i]))
        line += 1

cmap3d = cmap
cmap3d = None

# # crop above 1.2
# Xcross2 = []
# Ycross2 = []
# Zcross2 = []
# Xfan2 = []
# Yfan2 = []
# Zfan2 = []
# for i in range(len(Xcross)):
#     for j in range(len(Xcross[i])):
#         if Xcross[i, j] < 1.200001:
#             Xcross2[i].append(Xcross[i, j])
# print(Xcross.shape)
# print(Xcross)
Xcross = Xcross[:,:16]
Ycross = Ycross[:,:16]
Zcross = Zcross[:,:16]
Xfan = Xfan[:,:16]
Yfan = Yfan[:,:16]
Zfan = Zfan[:,:16]
# print(Xcross.shape)
# print(Xcross)
# print(Ycross)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf2 = ax.plot_surface(Xcross, Ycross, Zcross, cmap=cmap3d, color='m', #label='orthogonal',
                       linewidth=0, antialiased=False, vmin=0.0, vmax=0.05, alpha=0.6)
surf1 = ax.plot_surface(Xfan, Yfan, Zfan, cmap=cmap3d, #label='fan-like',
                       linewidth=0, antialiased=False, vmin=0.0, vmax=0.05, alpha=0.6)

ax.plot(mincross_x, betas, mincross, '-.', linewidth=2.0, color='red', label='orthogonal min')
ax.plot(minfan_x, betas, minfan, '-.', linewidth=2.0, color='black', label='fan-like min')
ax.plot(minint_x, betas, minint, '-.', linewidth=2.0, color='green', label='intersection')
plt.ylabel(r'springing angle ($\beta$)')
plt.xlabel(r'radius over length ($r/l_\mathrm{0}$)')
ax.set_xlim(0.5, 1.2)
# plt.zlabel(r'thickness over span ($t/s$)')
ax.set_zlabel(r'min thickness over span ($t_\mathrm{min}/s$)')
ax.legend()
ax.azim = 135
ax.dist = 9
ax.elev = 30
# fig.colorbar(surf1, shrink=0.5, aspect=5)
plt.show()


# # ---------------- By hand

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/pointed_crossvaut_intersection.csv'

intersec_min_results = [[], [], [], [], []]
intersec_min_x = [[], [], [], [], []]

with open(path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line = 0
    for row in csv_reader:
        if line > 1:
            for i in range(5):
                if row[1 + 2*i] == '':
                    continue
                intersec_min_results[i].append(float(row[1 + 2*i]))
                intersec_min_x[i].append(float(row[0 + 2*i]))
        line += 1

# print('\CROSS 3D')

degrees = []
i = 0
for serie in intersec_min_x:
    n = len(serie)
    degrees.append([i]*n)
    i += 10

Xint = np.array(intersec_min_x)
Yint = np.array(degrees)
Zint = np.array(intersec_min_results)

# print(Xint)
# print(Yint)
# print(Zint)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf1 = ax.plot_surface(Xint, Yint, Zint, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf1, shrink=0.5, aspect=5)
# plt.show()

# fig, ax = plt.subplots()
# CS = ax.contour(Xint, Yint, Zint, levels=20)
# plt.show()

# fig, ax = plt.subplots()
# CS = ax.contour(Yint, Xint, Zint, levels=20)
# plt.show()

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

# npts = 200

# Zint[4, 24] = 0.06
# Zint[0,0] = 0.0

xlim = (0.5, 1.2)
ylim = (0, 40)
ngridx = 200
ngridy = 100
# x = Xint.flatten()
# y = Yint.flatten()
# z = Zint.flatten()
x = Xfan.flatten()
y = Yfan.flatten()
z = Zfan.flatten()
# x = Xcross.flatten()
# y = Ycross.flatten()
# z = Zcross.flatten()
cmap = 'gist_rainbow'
# cmap = 'coolwarm_r'
# cmap = 'bwr'
cmap = 'seismic'
# cmap = 'BuGn_r'
inverse = True

# Values of CMAP: 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'

if inverse:
    x_ = y
    y = x
    x = x_
    xlim_ = ylim
    ylim = xlim
    xlim = xlim_
    ngridx_ = ngridy
    ngridy = ngridx
    ngridx_ = ngridx_

# x = np.random.uniform(-2, 2, npts)
# y = np.random.uniform(-2, 2, npts)
# z = x * np.exp(-x**2 - y**2)

npts = len(Xint.flatten())

print('test plot')
print(x.shape)
print(y.shape)
print(z.shape)
print(Xint.shape)
print(Xcross.shape)
print(Xfan.shape)
print(Zint.shape)
print(Zcross.shape)
print(Zfan.shape)

# --- -

vmin = 0.00
vmax = 0.05

print(Zint)

fig, ax2 = plt.subplots()

ax2.tricontour(x, y, z, levels=1000, linewidths=0.0, vmin=vmin, vmax=vmax)
cntr2 = ax2.tricontourf(x, y, z, levels=1000, cmap=cmap, vmin=vmin, vmax=vmax)

print(max(z), min(z))
# cntr2.set_clim(0.0, 0.05)
cbar = fig.colorbar(cntr2, ax=ax2, shrink=0.75, label=r'min thickness over span ($t_\mathrm{min}/s$)', ticks=[0.00 + 0.01*i for i in range(7)])  #, format=PercentFormatter(xmax=1, decimals=1))

# ax2.plot(x, y, 'ko', ms=3)
# ax2.set(xlim=(-2, 2), ylim=(-2, 2))
if inverse:
    ax2.plot(betas, mincross_x, '-.', linewidth=1.0, color='red', label='orthogonal min')
    ax2.plot(betas, minfan_x, '-.', linewidth=1.0, color='black', label='fan-like min')
    ax2.plot(betas, minint_x, '-.', linewidth=1.0, color='green', label='intersection')
else:
    ax2.plot(mincross_x, betas, 'o-')
    ax2.plot(minfan_x, betas, 'o-')
    ax2.plot(minint_x, betas, 'o-')
ax2.set(xlim=xlim, ylim=ylim)

ax2.set_aspect(40)
plt.xticks([0, 10, 20, 30, 40])
# plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.yticks([1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5])
plt.xlabel(r'springing angle ($\beta$)')
plt.ylabel(r'radius over length ($r/l_\mathrm{0}$)')
plt.legend()

plt.show()
