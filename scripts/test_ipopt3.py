import numpy as np
import scipy.sparse as sps
import ipopt
from scipy.optimize import approx_fprime


class Problem40(object):
    """ # Hock & Schittkowski test problem #40
            Basic structure  follows:
            - cyipopt example from https://pythonhosted.org/ipopt/tutorial.html#defining-the-problem
            - which follows ipopt's docs from: https://www.coin-or.org/Ipopt/documentation/node22.html
            Changes:
            - numerical-diff using scipy for function & constraints
            - removal of hessian-calculation
              - we will use limited-memory approximation
                - ipopt docs: https://www.coin-or.org/Ipopt/documentation/node31.html
              - (because i'm too lazy to reason about the math; lagrange and co.)
    """
    def __init__(self):
        self.num_diff_eps = 1e-8  # maybe tuning needed!

    def objective(self, x):
        # callback for objective
        return -np.prod(x)  # -x1 x2 x3 x4

    def constraint_0(self, x):
        return np.array([x[0]**3 + x[1]**2 -1])

    def constraint_1(self, x):
        return np.array([x[0]**2 * x[3] - x[2]])

    def constraint_2(self, x):
        return np.array([x[3]**2 - x[1]])

    def constraints(self, x):
        # callback for constraints
        return np.concatenate([self.constraint_0(x),
                               self.constraint_1(x),
                               self.constraint_2(x)])

    def gradient(self, x):
        # callback for gradient
        return approx_fprime(x, self.objective, self.num_diff_eps)

    def jacobian(self, x):
        # callback for jacobian
        return np.concatenate([
            approx_fprime(x, self.constraint_0, self.num_diff_eps),
            approx_fprime(x, self.constraint_1, self.num_diff_eps),
            approx_fprime(x, self.constraint_2, self.num_diff_eps)])

    def hessian(self, x, lagrange, obj_factor):
        return False  # we will use quasi-newton approaches to use hessian-info

    # progress callback
    def intermediate(
            self,
            alg_mod,
            iter_count,
            obj_value,
            inf_pr,
            inf_du,
            mu,
            d_norm,
            regularization_size,
            alpha_du,
            alpha_pr,
            ls_trials
            ):

        print("Objective value at iteration #%d is - %g" % (iter_count, obj_value))

# Remaining problem definition; still following official source:
# http://www.ai7.uni-bayreuth.de/test_problem_coll.pdf

# start-point -> infeasible
x0 = [0.8, 0.8, 0.8, 0.8]

# variable-bounds -> empty => np.inf-approach deviates from cyipopt docs!
lb = [-np.inf, -np.inf, -np.inf, -np.inf]
ub = [np.inf, np.inf, np.inf, np.inf]

# constraint bounds -> c == 0 needed -> both bounds = 0
cl = [0, 0, 0]
cu = [0, 0, 0]

p40 = Problem40()
jac = p40.jacobian(x0)
grad = p40.gradient(x0)
print('Gradient')
print(grad)
print(grad.shape)
print('Jacobian')
print(jac)
print(jac.shape)

nlp = ipopt.problem(
            n=len(x0),
            m=len(cl),
            problem_obj=Problem40(),
            lb=lb,
            ub=ub,
            cl=cl,
            cu=cu
            )

# IMPORTANT: need to use limited-memory / lbfgs here as we didn't give a valid hessian-callback
nlp.addOption(b'hessian_approximation', b'limited-memory')
x, info = nlp.solve(x0)
print(x)
print(info)

# CORRECT RESULT & SUCCESSFUL STATE
