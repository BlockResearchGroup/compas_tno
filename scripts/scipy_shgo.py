from scipy.optimize import shgo
import numpy as np

def f(x):  # (cattle-feed)
    return 24.55*x[0] + 26.75*x[1] + 39*x[2] + 40.50*x[3]

def g1(x):
    return 2.3*x[0] + 5.6*x[1] + 11.1*x[2] + 1.3*x[3] - 5  # >=0

def g2(x):
    return (12*x[0] + 11.9*x[1] +41.8*x[2] + 52.1*x[3] - 21
            - 1.645 * np.sqrt(0.28*x[0]**2 + 0.19*x[1]**2
                            + 20.5*x[2]**2 + 0.62*x[3]**2))# >=0

def h1(x):
    return x[0] + x[1] + x[2] + x[3] - 1  # == 0

cons = [{'type': 'ineq', 'fun': g1},
        {'type': 'ineq', 'fun': g2},
        {'type': 'eq', 'fun': h1}]
bounds = [(0, 1.0),]*4
res = shgo(f, bounds, n=3, iters=3, constraints=cons, options={'disp': True})
print(res)
x = res['x']

print(g1(x))
