import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
Y = np.arange(0.5, 1.00001, 0.05)
X = np.arange(0, 41, 10)
print(len(X))
print(len(Y))
X, Y = np.meshgrid(X, Y)
# print(Y.shape)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

# cross
Zcross = [0.335,	0.303,	0.231,	0.154,	0.096,
0.301,	0.261,	0.177,	0.095,	0.075,
0.275,	0.232,	0.145,	0.139,	0.154,
0.259,	0.214,	0.185,	0.200,	0.233,
0.249,	0.229,	0.239,	0.259,	0.304,
0.274,	0.275,	0.287,	0.321,	0.370,
0.316,	0.318,	0.331,	0.378,	0.429,
0.354,	0.357,	0.372,	0.431,	0.485,
0.389,	0.393,	0.409,	0.480,	0.537,
0.422,	0.426,	0.445,	0.527,	0.587,
0.452,	0.457,	0.478,	0.571,	0.634]

# fan
Zfan = [0.469,	0.442,	0.374,	0.301,	0.233,
0.433,	0.402,	0.326,	0.241,	0.165,
0.413,	0.379,	0.295,	0.204,	0.136,
0.397,	0.359,	0.269,	0.173,	0.138,
0.384,	0.343,	0.245,	0.168,	0.189,
0.374,	0.328,	0.228,	0.209,	0.235,
0.364,	0.315,	0.227,	0.249,	0.277,
0.357,	0.307,	0.259,	0.287,	0.320,
0.352,	0.299,	0.288,	0.321,	0.360,
0.348,	0.299,	0.318,	0.354,	0.398,
0.345,	0.323,	0.346,	0.385,	0.434]

Zmixed = [0.411, 0.385,	0.318,	0.241,	0.173,
0.377,	0.347,	0.285,	0.189,	0.115,
0.359,	0.321,	0.236,	0.168,	0.193,
0.345,	0.305,	0.224,	0.231,	0.251,
0.337,	0.294,	0.270,	0.292,	0.315,
0.332,	0.295,	0.307,	0.348,	0.375,
0.331,	0.336,	0.354,	0.397,	0.428,
0.368,	0.374,	0.398,	0.434,	0.474,
0.401,	0.408,	0.438,	0.478,	0.524,
0.432,	0.440,	0.476,	0.520,	0.579,
0.461,	0.470,	0.512,	0.562,	0.623]

z_mixed = []

for i in range(len(Zcross)):
    z_mixed.append(min(Zcross[i], Zfan[i]))

# Z = np.array(Zcross).reshape(11, 5)
# Z = np.array(Zfan).reshape(11, 5)
Z = np.array(Zmixed).reshape(11, 5)
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
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, levels=20)
plt.show()
