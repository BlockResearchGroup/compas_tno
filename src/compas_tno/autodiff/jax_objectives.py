from jax.numpy import dot
from jax.numpy import diag
from jax.numpy import zeros
from jax.numpy import array
from jax.numpy.linalg import norm
from jax import grad

from compas_tno.algorithms import q_from_variables

from compas_tno.autodiff.jax_equilibrium import xyz_from_q_jax


def objective_selector_jax(objective):
    """Select objective callable and gradient vector based on the desired objective function using JAX

    Parameters
    ----------
    objective : str
        The name of the objective. See ``Optimiser`` for a complete list of objectivees.

    Returns
    -------
    fobj : callable
        Callable to compute the value of the objective function in the point.
    fgrad : callable
        Callable to compute the gradient of the objective function in the point.

    """

    if objective == 'min':
        fobj = f_min_thrust_jax
        fgrad = grad(f_min_thrust_jax)
    elif objective == 'max':
        fobj = f_max_thrust_jax
        fgrad = grad(f_max_thrust_jax)
    elif objective == 'loadpath':
        # fgrad = grad(f_loadpath_general)
        raise NotImplementedError
    elif objective == 'target' or objective == 'bestfit':
        # fgrad = grad(f_bestfit)
        raise NotImplementedError
    elif objective == 'feasibility':
        # fgrad = grad(f_constant)
        raise NotImplementedError
    elif objective == 'hor_projection':
        # fgrad = grad(f_horprojection)
        raise NotImplementedError
    elif objective == 't':  # analytical reduce thickness
        # fgrad = grad(f_reduce_thk)
        raise NotImplementedError
    elif objective == 's':  # tight UB and LB 0 -> 1/2
        # fgrad = grad(f_tight_crosssection)
        raise NotImplementedError
    elif objective == 'n':  # vector n offset the surfaces -> larger the better (higher GSF)
        # fgrad = grad(f_tight_crosssection)
        raise NotImplementedError
    elif objective == 'lambdh':  # vector lambda as hor multiplier larger the better (higher GSF)
        # fgrad = grad(f_tight_crosssection)
        raise NotImplementedError
    elif objective == 'Ecomp-linear':  # vector lambda as hor multiplier larger the better (higher GSF)
        # fgrad = grad(f_complementary_energy)
        raise NotImplementedError
    elif objective == 'Ecomp-nonlinear':
        # fgrad = grad(f_complementary_energy_nonlinear)
        raise NotImplementedError
    elif objective == 'max_section':
        # fgrad = grad(f_max_section)
        raise NotImplementedError
    else:
        print('Please, provide a valid objective for the optimisation')
        raise NotImplementedError

    return fobj, fgrad


def f_min_thrust_jax(variables, M):
    """Objective function to minimise the horizontal thrust ready for JAX

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    update_geometry = False  # In the typical case the geometry does not need to be updated
    k = M.k
    nb = len(M.fixed)

    P_Xh_fixed = M.P[M.fixed][:, :2]  # Horizontal loads in the fixed vertices
    P_free = M.P[M.free]  # Loads in the free vertices

    X = array(M.X)
    Xh = X[:, :2]  # xy of the thrust (i.e. horizontal projection)

    SPLU_D = None

    qid = variables[:k].reshape(-1, 1)
    q = q_from_variables(qid, M.B, M.d)
    Q = diag(q.flatten())

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb].reshape(-1, 2, order='F')
        X.at[M.fixed, :2].set(xyb)
        update_geometry = True
    if 'zb' in M.variables:
        zb = variables[-nb:]
        X.at[M.fixed, 2].set(zb.flatten())
        if 'fixed' not in M.features:
            update_geometry = True

    if update_geometry:
        if 'update-loads' in M.features:
            pass

        X[M.free] = xyz_from_q_jax(q, P_free, X[M.fixed], M.Ci, M.Cit, M.Cb, SPLU_D=SPLU_D)
        Xh = X[:, :2]

    CfQC = dot(dot(M.Cb.transpose().toarray(), Q), M.C.toarray())
    Rh = dot(CfQC, Xh) - P_Xh_fixed
    f = sum(norm(Rh, axis=1))

    return f


def f_max_thrust_jax(variables, M):
    """Objective function to maximise the horizontal thrust ready for JAX.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    return -1 * f_min_thrust_jax(variables, M)
