from jax.numpy import dot
from jax.numpy import diag
from jax.numpy import array
from jax.numpy import vstack
from jax.numpy import hstack

from jax.scipy.linalg import lu_solve
from jax.scipy.linalg import lu_factor

from compas_tno.autodiff.jax_swt import weights_from_xyz_jax


def xyz_from_q_jax(q, Pi, Xb, Ci, Cit, Cb, update_loads=False, M=None):
    """ Calculate coordinates's xyz the q's.

    Parameters
    ----------
    q : array [m x 1]
        Force densities of all edges.
    Pi : array [ni x 3]
        External forces applied to the free vertices.
    Xb : array [nb x 3]
        Position of the fixed vertices.
    Ci : array [m x ni]
        Connectivity matrix on the free vertices
    Cit : array [ni x m]
        Transpose of connectivity matrix on the free vertices
    Cb : array [m x nb]
        Connectivity matrix on the fixed vertices
    update_loads : bool, optional
        Whether or not vertical loads are updated according to the vertex tributary area, by default False

    Returns
    -------
    Xfree : array [ni x 3]
        x, y, z coordinates of the nodes

    """
    Cit = Cit.toarray()
    Ci = Ci.toarray()
    Cb = Cb.toarray()
    Q = diag(q.flatten())

    CiQCb = dot(dot(Cit, Q), Cb)  # Cit.dot(diags(q.flatten())).dot(Cb)
    CiQCi = dot(dot(Cit, Q), Ci)  # Cit.dot(diags(q.flatten())).dot(Ci)
    b = Pi - dot(CiQCb, Xb)  # Pi - CiQCb.dot(Xb)
    lu_fac = lu_factor(CiQCi)
    Xfree = lu_solve(lu_fac, b)

    if update_loads:
        if M is None:
            raise ValueError('Please provide Problem object (M)')
        X = dot(M.Pmatrix, vstack(Xfree, Xb))
        pz = weights_from_xyz_jax(X, M.F, M.V0, M.V1, M.V2)
        Pi = dot(M.Pmatrix[:, :len(M.free)], pz[array(M.free)])
        b = Pi - dot(CiQCb, Xb)
        Xfree = lu_solve(lu_fac, b)

    return Xfree


def z_from_q_jax(q, pz, zb, Ci, Cit, Cb, update_loads=False, M=None):
    """ Calculate coordinates's xyz the q's.

    Parameters
    ----------
    q : array [m x 1]
        Force densities of all edges.
    pz : array [ni x 1]
        External forces applied to the free vertices.
    zb : array [nb x 1]
        Position of the fixed vertices.
    Ci : array [m x ni]
        Connectivity matrix on the free vertices
    Cit : array [ni x m]
        Transpose of connectivity matrix on the free vertices
    Cb : array [m x nb]
        Connectivity matrix on the fixed vertices
    update_loads : bool, optional
        Whether or not vertical loads are updated according to the vertex tributary area, by default False

    Returns
    -------
    zfree : array [ni x 1]
        z coordinates of the free nodes

    """
    Cit = Cit.toarray()
    Ci = Ci.toarray()
    Cb = Cb.toarray()
    Q = diag(q.flatten())

    CiQCb = dot(dot(Cit, Q), Cb)  # Cit.dot(diags(q.flatten())).dot(Cb)
    CiQCi = dot(dot(Cit, Q), Ci)  # Cit.dot(diags(q.flatten())).dot(Ci)
    b = pz - dot(CiQCb, zb)  # Pi - CiQCb.dot(Xb)
    lu_fac = lu_factor(CiQCi)
    zfree = lu_solve(lu_fac, b)

    if update_loads:
        if M is None:
            raise ValueError('Please provide Problem object (M)')
        Xfixed = M.X[M.fixed]
        Xb = hstack([Xfixed[:2], zb])
        X = dot(M.Pmatrix, vstack(Xfree, Xb))
        pz = weights_from_xyz_jax(X, M.F, M.V0, M.V1, M.V2)
        Pi = dot(M.Pmatrix[:, :len(M.free)], pz[array(M.free)])
        b = Pi - dot(CiQCb, Xb)
        Xfree = lu_solve(lu_fac, b)

    return zfree


def q_from_variables_jax(qid, B, d, lambd=1.0):
    r""" Calculate q's from the force parameters (independent edges) using JAX.

    Parameters
    ----------
    qid : array [kx1]
        Force density vector.
    B : array [m x k]
        The linear map on the force densities.
    d : array [mx1]
        A particular solution of the equilibrium.
    lambd : float, optional
        Lambda multiplier applied to horizontal loads. The default is 1.0.

    Returns
    -------
    q : array [m x 1]
        Force densities on all edges.

    Notes
    -------
    This function works for load multiplers and is prefered. To see details about the implementation check ``q_from_qid``.

    Reference
    ---------
    Block and Lachauer, 2014...

    """

    q = dot(B, qid) + lambd * d

    return q
