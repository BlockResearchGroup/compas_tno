from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_independents
from compas_tno.algorithms import initialise_problem
from compas_tno.algorithms import z_from_form
import compas_tno

from compas.numerical import normrow

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags

import torch as th
from torch.autograd.gradcheck import zero_gradients

from numpy import zeros
from numpy import identity
from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide


# -------------- NUMPY CALCULATIONS


def q_from_qid_trad(qid, Ei, Edinv, p, ind, dep, m):
    qdep = Edinv.dot(- p + Ei.dot(qid))
    q = zeros((m, 1))
    q[dep] = qdep
    q[ind] = qid
    return q


def z_from_qid_trad(qid, Ei, Edinv, p, ind, dep, m, Ci, Cit, zfixed, pzfree):
    q = q_from_qid_trad(qid, Ei, Edinv, p, ind, dep, m)
    zi = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pzfree - Cit.dot(diags(q.flatten())).dot(Cf).dot(zfixed))
    return zi


def reac_bound_trad(qid, Ei, Edinv, p, ind, dep, m, C, Cf, zfixed, pfixed, xyz):
    q = q_from_qid_trad(qid, Ei, Edinv, p, ind, dep, m)
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    R = CfQC.dot(xyz) - pfixed
    length = abs(multiply(zfixed, divide(R[:, 0], R[:, 2]).reshape(-1, 1)))
    return length


# -------------- SIMPLIFIED MATRICES


def q_from_qid_matrix(qid, p, A_Ik, B_Ik):
    q = A_Ik.dot(p) + B_Ik.dot(qid)
    return q


def z_from_qid_matrix(qid, p, A_Ik, B_Ik, Ci, Cit, zfixed):
    q = A_Ik.dot(p) + B_Ik.dot(qid)
    zi = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(zfixed))
    return zi


def reac_bound_matrix(qid, p, A_Ik, B_Ik, C, Cf, zfixed, pfixed, xyz):
    q = q_from_qid_matrix(qid, p, A_Ik, B_Ik)
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    R = CfQC.dot(xyz) - pfixed
    length = abs(multiply(zfixed, divide(R[:, 0], R[:, 2]).reshape(-1, 1)))
    return length


def f_min_thrust_matrix(qid, p, A_Ik, B_Ik, C, Cf, xy):
    q = q_from_qid_matrix(qid, p, A_Ik, B_Ik)
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    Rh = CfQC.dot(xy)
    f = sum(normrow(Rh))
    return f

# -------------- PYTORCH IMPLEMENTATION


def q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik):
    q = th.mm(A_Ik, p) + th.mm(B_Ik, qid)
    return q


def z_from_qid_matrix_pytorch(qid, p, A_Ik, B_Ik, Ci, Cit, Cf, zfixed, pzfree):
    q = th.mm(A_Ik, p) + th.mm(B_Ik, qid)
    Q = th.diagflat(q)
    Ai = th.mm(th.mm(Cit, Q), Ci)
    Af = th.mm(th.mm(Cit, Q), Cf)
    b = pzfree - th.mm(Af, zfixed)
    X, LU = th.solve(b, Ai)
    zi = X
    return zi


def reac_bound_matrix_pytorch(qid, p, A_Ik, B_Ik, C, Cf, zfixed, pfixed, xyz):
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = th.diagflat(q)
    CfQC = th.mm(th.mm(Cf.t(), Q), C)
    R = th.mm(CfQC, xyz) - pfixed
    length_x = abs(th.mul(zfixed, th.div(R[:, 0], R[:, 2]).reshape(-1, 1)))
    length_y = abs(th.mul(zfixed, th.div(R[:, 1], R[:, 2]).reshape(-1, 1)))
    length = th.cat([length_x, length_y])
    return length

# -------------- JACOBIAN

def compute_jacobian(inputs, outputs):
    d_otp = outputs.size()[0]
    d_inp = inputs.size()[0]
    jacobian = th.zeros(d_otp, d_inp)
    grad_output = th.zeros(d_otp, 1)
    for i in range(d_otp):
        zero_gradients(inputs)
        grad_output.zero_()
        grad_output[i] = 1
        outputs.backward(grad_output, retain_graph=True)
        jacobian[i] = inputs.grad.data.flatten()
    return th.transpose(jacobian, dim0=0, dim1=1)


# -------------- OBJECTIVE FUNCTIONS

def f_min_thrust_pytorch(qid, p, A_Ik, B_Ik, C, Cf, xy):
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = th.diagflat(q)
    CfQC = th.mm(th.mm(Cf.t(), Q), C)
    Rh = th.mm(CfQC, xy)
    R = th.norm(Rh, dim=1)
    f = th.sum(R)
    return f


def f_max_thrust_pytorch(qid, p, A_Ik, B_Ik, C, Cf, xy):
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = th.diagflat(q)
    CfQC = th.mm(th.mm(Cf.t(), Q), C)
    Rh = th.mm(CfQC, xy)
    R = th.norm(Rh, dim=1)
    f = -1 * th.sum(R)
    return f


# -------------- GRADIENT OF OBJECTIVE FUNCTIONS


def compute_grad(variables, f):
    f.backward(retain_graph=True)
    grad = variables.grad.data
    return grad

# ==============================================================================
# Main
# ==============================================================================


if __name__ == "__main__":

    file_address = compas_tno.get('test.json')
    form = FormDiagram.from_json(file_address)
    form = z_from_form(form)
    # plot_independents(form).show()

    # Calculate independents and get back this zillions of matrices
    args = initialise_problem(form)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    # Prepare the matrix to obtain a simple function from R(k) -> R(m)
    # That get the independents (size k) and return the complete vector of q (size m)

    m = len(q)
    B_Ik = zeros((m, len(ind)))
    A_Ik = zeros((m, len(p)))
    B_Ik[ind, :] = identity(len(ind))
    B_Ik[dep, :] = Edinv*Ei
    A_Ik[dep, :] = -1*Edinv.toarray()

    # Create the tensors to use pytorch

    A_Ik_th = th.tensor(A_Ik)
    B_Ik_th = th.tensor(B_Ik)
    p_th = th.tensor(p)

    # Define my variables - the ones I care to derivate with regards to...

    qid_th = th.tensor(q[ind], requires_grad=True)
    zb_th = th.tensor(z[fixed], requires_grad=True)

    # Compute the q's as tensor

    q_ = q_from_qid_matrix_pythorch(qid_th, p_th, A_Ik_th, B_Ik_th)

    # Compute the zi's as a tensor

    C_th = th.tensor(C.toarray())
    Ci_th = th.tensor(Ci.toarray())
    Cit_th = Ci_th.t()
    Cf_th = th.tensor(Cf.toarray())
    pzfree = th.tensor(pz[free])

    zi_th = z_from_qid_matrix_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, Ci_th, Cit_th, Cf_th, zb_th, pzfree)

    # Compute Reaction Bounds

    xyz = th.tensor(hstack([x, y, z]))
    xy = th.tensor(hstack([x, y]))
    pfixed = th.tensor(hstack([px, py, pz])[fixed])

    length = reac_bound_matrix_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, C_th, Cf_th, zb_th, pfixed, xyz)

    # Compute Jacobian from q's based on qid by mark

    jac_q_qid = compute_jacobian(qid_th, q_)
    print(len(q_))
    print(jac_q_qid.shape)

    # Compute Jacobian from q's based on zb by mark (This has to be all zero)

    # jac_q_zb = compute_jacobian(zb_th, q_)
    # print(len(q_))
    # print(jac_q_zb.shape)

    # Compute Jacobian from zi's based on qid and zb by mark

    jac_z_qid = compute_jacobian(qid_th, zi_th)
    print(len(zi_th))
    print(jac_z_qid.shape)

    jac_z_zb = compute_jacobian(zb_th, zi_th)
    print(len(zi_th))
    print(jac_z_zb.shape)

    # Compute Jacobian from reac_bounds based on qid and zb by mark

    jac_length_qid = compute_jacobian(qid_th, length)
    print(len(length))
    print(jac_length_qid.shape)

    jac_length_zb = compute_jacobian(zb_th, length)
    print(len(length))
    print(jac_length_zb.shape)

    # Compute Objective Function and Grad (MIN)

    f = f_min_thrust_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, C_th, Cf_th, xy)
    print('f = ', f)

    grad_f = compute_grad(qid_th, f)
    print(grad_f)

    # Compute Objective Function and Grad (MAX)

    f = f_max_thrust_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, C_th, Cf_th, xy)
    print('f = ', f)

    grad_f = compute_grad(qid_th, f)
    print(grad_f)

    # Compute Objective Function and Grad (LOADPATH)

    # f = f_bestfit_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, C_th, Cf_th, xy)
    # print('f = ', f)

    # grad_f = compute_grad(qid_th, f)
    # print(grad_f)

    # Compute Objective Function and Grad (BESTFIT)

    # f = f_bestfit_pytorch(qid_th, p_th, A_Ik_th, B_Ik_th, C_th, Cf_th, xy)
    # print('f = ', f)

    # grad_f = compute_grad(qid_th, f)
    # print(grad_f)



    # Set Optimisation



