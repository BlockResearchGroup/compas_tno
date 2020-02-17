import numpy as np
import matplotlib.pyplot as plt
import math

density = 200
number_lines = 10
rmax = 18

delta = 2*math.pi/density
theta = np.arange(0, 3/2 * math.pi + delta, delta)
ks = np.arange(0, number_lines, 1)

ax = plt.subplot(111, projection='polar')

# Draw segments of compression lines

for k in ks:
    theta_n = [[]]
    r_n = [[]]
    for i in range(1, len(theta)):
        n2 = 2/3*k/math.sin(2/3*theta[i])
        if n2 >= 0:
            n = math.pow(n2, 3/2)
        else:
            n = math.pow(-1*n2, 3/2)
        plot = False
        if 0 < theta[i] < 3/2 * math.pi:
            quad = 0
            plot = True
        if plot == True:
            theta_n[quad].append(theta[i])
            r_n[quad].append(n)
    for i in range(len(theta_n)):
        ax.plot(theta_n[i], r_n[i], color='red')

# Draw segments of tension lines
# theta = np.roll(theta, -1 * int(len(theta)/8+1))

for k in ks:
    theta_m = [[],[]]
    r_m = [[],[]]
    for i in range(len(theta)):
        m2 = 2/3*k/math.cos(2/3*theta[i])
        if m2 >= 0:
            m = math.pow(m2, 3/2)
        else:
            m = math.pow(-1*m2, 3/2)
        plot = False
        if 0 <= theta[i] < 3/8 * 2*math.pi:
            quad = 0
            plot = True
        elif 3/8 * 2*math.pi < theta[i] <= 3/2 * math.pi:
            quad = 1
            plot = True
        if plot == True:
            theta_m[quad].append(theta[i])
            r_m[quad].append(m)
    for i in range(len(theta_m)):
        ax.plot(theta_m[i], r_m[i], color='blue')

ax.set_rmax(rmax)
ax.set_rticks([0.0])  # Less radial ticks
# ax.set_thetaticks([0.0])
# ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)
ax.set_title("Maxwell Force Diagram", va='bottom')
plt.show()
