import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(0.5, 1.00001, 0.05)
Y = np.arange(0, 41, 5)
print(len(X))
print(len(Y))
X, Y = np.meshgrid(X, Y)
# print(Y.shape)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

Xcross = [[0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
          [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]
          ]

Ycross = [[0]*11, [5]*11, [10]*11, [15]*11, [20]*11, [25]*11, [30]*11, [35]*11, [40]*11]

# cross
Zcross = [[
3.34671137515692,2.98093292008925,2.73957895366525,2.57714490820888,2.46466874015549,2.50616540391259,2.57465643823119,2.73140657040956,2.90970390993919,3.07071168747796,3.19257600528445,
3.2673132814599,2.8866803452628,2.62999826751377,2.46040587970138,2.33495256078826,2.3771904050686,2.46127563604126,2.65377007725544,2.84406526075481,3.05698334109704,3.15305833514612,
3.16999251318562,2.60880522633673,2.31336547018319,2.12781331108047,2.08919776736852,2.17052068545024,2.38842889745827,2.61479987059196,2.79606433579081,2.83921486133035,2.89859761360004,
2.68445754698909,2.21457803402135,1.89373179575741,1.75926490464664,1.83057440921253,2.10041299942886,2.36051546246387,2.48248997657753,2.56903067240417,2.77762572198887,2.96575662485103,
2.28994661819171,1.77263959259923,1.44324687003251,1.48576765130317,1.78022715409977,2.05296935438907,2.13903746178817,2.405637018108,2.64336024522949,2.87730648796638,3.09542692753,
1.89017121150722,1.33370563498001,1.12542230681157,1.42700106049621,1.74304201018361,1.9057600789049,2.21263510428331,2.51277177870023,2.78769915822948,3.03570017338917,3.26079363526039,
1.51406724185019,0.951331779943179,1.03708607960665,1.40841645565269,1.62385759137973,1.9883765269099,2.34114266484981,2.65499499425379,2.93639362288847,3.19041858314757,3.4211225201339,
1.17765959158888,0.678424658287976,1.02281425888703,1.32217527733985,1.70530588102899,2.11208444416089,2.47100325593845,2.79046892372172,3.07878030476081,3.36592933327307,3.62706733121254,
0.888074133291907,0.570522025158226,0.977755714183095,1.34442094519572,1.81623936776764,2.22847979959217,2.5965776004978,2.95305839382457,3.27308742670292,3.56221815232204,3.82490995829504,
]]


# # fan
# Zfan = [0.469,	0.442,	0.374,	0.301,	0.233,
#         0.433,	0.402,	0.326,	0.241,	0.165,
#         0.413,	0.379,	0.295,	0.204,	0.136,
#         0.397,	0.359,	0.269,	0.173,	0.138,
#         0.384,	0.343,	0.245,	0.168,	0.189,
#         0.374,	0.328,	0.228,	0.209,	0.235,
#         0.364,	0.315,	0.227,	0.249,	0.277,
#         0.357,	0.307,	0.259,	0.287,	0.320,
#         0.352,	0.299,	0.288,	0.321,	0.360,
#         0.348,	0.299,	0.318,	0.354,	0.398,
#         0.345,	0.323,	0.346,	0.385,	0.434]

# Zmixed = [0.411, 0.385,	0.318,	0.241,	0.173,
#         0.377,	0.347,	0.285,	0.189,	0.115,
#         0.359,	0.321,	0.236,	0.168,	0.193,
#         0.345,	0.305,	0.224,	0.231,	0.251,
#         0.337,	0.294,	0.270,	0.292,	0.315,
#         0.332,	0.295,	0.307,	0.348,	0.375,
#         0.331,	0.336,	0.354,	0.397,	0.428,
#         0.368,	0.374,	0.398,	0.434,	0.474,
#         0.401,	0.408,	0.438,	0.478,	0.524,
#         0.432,	0.440,	0.476,	0.520,	0.579,
#         0.461,	0.470,	0.512,	0.562,	0.623]


# Find the minimum among them

# z_mixed = []

# for i in range(len(Zcross)):
#     z_mixed.append(min(Zcross[i], Zfan[i]))

# Z = np.array(Zcross).reshape(11, 5)
Z = np.array(Zcross).reshape(9, 11)
# Z = np.array(Zfan).reshape(11, 5)
# Z = np.array(Zmixed).reshape(11, 5)
# Z = np.array(z_mixed).reshape(11, 5)

print(X.shape)
print(Y.shape)
print(Z.shape)
print(X)
print(Y)
print(Z)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

# Basic contour plot
fig, ax = plt.subplots(figsize=(12, 5))
# im = ax.imshow(Z, interpolation='bilinear', origin='lower',
#                cmap=cm.gray)
# im = ax.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
            #    origin='lower', extent=[0.5, 1.0, 0, 40])
            #    vmax=abs(Z).max(), vmin=-abs(Z).max())

levels = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

CS = ax.contour(X, Y, Z, levels=levels, colors='k')
ax.clabel(CS, levels[2::2], inline=True, fontsize=10)

CS = ax.contourf(X, Y, Z, cmap='seismic', levels=1000, vmin=0.00, vmax=5.0)
cbar = fig.colorbar(CS, ax=ax, shrink=0.75, label=r'($t_\mathrm{min}/s$)', ticks=[0.00 + 0.5*i for i in range(8)])

# CBI = fig.colorbar(im, orientation='horizontal', shrink=0.8)
plt.show()


### TRY the colors

vmin = 0.00
vmax = 0.05
cmap = 'seismic'

fig, ax2 = plt.subplots()

X = X.reshape(-1, 1)
Y = Y.reshape(-1, 1)
Z = Z.reshape(-1, 1)

ax2.tricontour(X, Y, Z, levels=1000, linewidths=0.0, vmin=vmin, vmax=vmax)
cntr2 = ax2.tricontourf(X, Y, Z, levels=1000, cmap=cmap, vmin=vmin, vmax=vmax)

print(max(z), min(z))
# cntr2.set_clim(0.0, 0.05)
cbar = fig.colorbar(cntr2, ax=ax2, shrink=0.75, label=r'min thickness over span ($t_\mathrm{min}/s$)', ticks=[0.00 + 0.01*i for i in range(7)])  #, format=PercentFormatter(xmax=1, decimals=1))

# ax2.plot(x, y, 'ko', ms=3)
# ax2.set(xlim=(-2, 2), ylim=(-2, 2))
# if inverse:
#     ax2.plot(betas, mincross_x, '-.', linewidth=1.0, color='red', label='orthogonal min')
#     ax2.plot(betas, minfan_x, '-.', linewidth=1.0, color='black', label='fan-like min')
#     ax2.plot(betas, minint_x, '-.', linewidth=1.0, color='green', label='intersection')
# else:
#     ax2.plot(mincross_x, betas, 'o-')
#     ax2.plot(minfan_x, betas, 'o-')
#     ax2.plot(minint_x, betas, 'o-')
# ax2.set(xlim=xlim, ylim=ylim)

ax2.set_aspect(40)
plt.xticks([0, 10, 20, 30, 40])
# plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.yticks([1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5])
plt.xlabel(r'springing angle ($\beta$)')
plt.ylabel(r'radius over length ($r/l_\mathrm{0}$)')
plt.legend()

plt.show()
