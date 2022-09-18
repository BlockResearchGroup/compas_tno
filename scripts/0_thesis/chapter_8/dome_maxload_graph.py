import matplotlib.pyplot as plt
import numpy as np

def box_limits(box):
    return [box.x0, box.y0, box.width*0.8, box.height]

size_plots = (10, 5)

size_axis_label = 12
size_axis_data = 12
size_legend = 12

#### MAX LOAD APEX

# Plot for graph of curves and discretisations
series = []
x = []
y = []
styles = []
colours = []

radius = 5.0
bench = 0.042
discr = [4, 8, 12, 16, 20, 24]

# meridians = 12
series.append(r'$n_\mathrm{M} = 12$')
y.append(np.array([0.1773, 0.1529, 0.1468, 0.1443, 0.1430, 0.1421]))

# meridians = 16
series.append(r'$n_\mathrm{M} = 16$')
y.append(np.array([0.1773, 0.1527, 0.1465, 0.1439, 0.1426, 0.1419]))

# meridians = 20
series.append(r'$n_\mathrm{M} = 20$')
y.append(np.array([0.1769, 0.1526, 0.1464, 0.1438, 0.1425, 0.1417]))

# meridians = 24
series.append(r'$n_\mathrm{M} = 24$')
y.append(np.array([0.1768, 0.1525, 0.1463, 0.1437, 0.1424, 0.1416]))

# for i in range(len(series)):
#     plt.plot(x[i], y[i], 'o-', label=series[i])


fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
ax = plt.axes()
ax.plot(discr, y[0], 'o-', label=series[0])
ax.plot(discr, y[1], 'o-', label=series[1])
ax.plot(discr, y[2], 'o-', label=series[2])
ax.plot(discr, y[3], 'o-', label=series[3])

ax.legend(fontsize=size_legend)
ax.set_xlabel(r'$n_\mathrm{P}$', size=size_axis_label)
ax.set_ylabel(r'$P_\mathrm{max}/W$', size=size_axis_label)
ax.tick_params(labelsize=size_axis_data)
ax.set_xlim(4-1, 24+1)
ax.set_ylim(0.13, 0.18)
ax.set_xticks([4, 8, 12, 16, 20, 24])

box = ax.get_position()
ax.set_position(box_limits(box))

plt.show()

#### MAX LOAD OFF-CENTERED

size_plots = (10, 5)

size_axis_label = 12
size_axis_data = 12
size_legend = 12

# Plot for graph of curves and discretisations
series = []
x = []
y = []
styles = []
colours = []

radius = 5.0
bench = 0.042
discr = [4, 8, 12, 16, 20, 24]

# meridians = 12
series.append(r'$n_\mathrm{M} = 12$')
y.append(np.array([0.0752892152655848, 0.068904147238022, 0.0628975202845967, 0.0602227346228325, 0.0597086182266271, 0.059551755421957]))

# meridians = 16
series.append(r'$n_\mathrm{M} = 16$')
y.append(np.array([0.0687844860251164, 0.0604887590199928, 0.0551425335981768, 0.0531115774822682, 0.0526420109886243, 0.0524245173469545]))

# meridians = 20
series.append(r'$n_\mathrm{M} = 20$')
y.append(np.array([0.063378259728144, 0.0559250599471017, 0.0513215438782695, 0.0493721421814023, 0.0489747283086051, 0.0487948356968229]))

# meridians = 24
series.append(r'$n_\mathrm{M} = 24$')
y.append(np.array([0.0600988928861676, 0.0529583846714335, 0.0491775022260183, 0.0474341922360024, 0.0470010443135512, 0.0467568279195238]))

# for i in range(len(series)):
#     plt.plot(x[i], y[i], 'o-', label=series[i])


fig = plt.figure(figsize=size_plots)  # try 12, 4
# ax = plt.subplot(111)
ax = plt.axes()
ax.plot(discr, y[0], 'o-', label=series[0])
ax.plot(discr, y[1], 'o-', label=series[1])
ax.plot(discr, y[2], 'o-', label=series[2])
ax.plot(discr, y[3], 'o-', label=series[3])

ax.legend(fontsize=size_legend)
ax.set_xlabel(r'$n_\mathrm{P}$', size=size_axis_label)
ax.set_ylabel(r'$P_\mathrm{max}/W$', size=size_axis_label)
ax.tick_params(labelsize=size_axis_data)
ax.set_xlim(4-1, 24+1)
ax.set_ylim(0.03, 0.08)
ax.set_xticks([4, 8, 12, 16, 20, 24])

box = ax.get_position()
ax.set_position(box_limits(box))

plt.show()
