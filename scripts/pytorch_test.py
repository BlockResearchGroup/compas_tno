from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_independents
from compas_tno.algorithms import initialise_problem
import compas_tno

import torch as th
from torch.autograd.gradcheck import zero_gradients

from numpy import zeros
from numpy import identity

def q_from_qid_prod(qid, Ei, Edinv, p, ind, dep, m):
    B_Ik = zeros((m, len(ind)))
    A_Ik = zeros((m, len(p)))
    B_Ik[ind, :] = identity(len(ind))
    B_Ik[dep, :] = Edinv*Ei
    A_Ik[dep, :] = -1*Edinv.toarray()
    q = A_Ik.dot(p) + B_Ik.dot(qid)
    return q

def q_from_qid_matrix(qid, p, A_Ik, B_Ik):
    q = A_Ik.dot(p) + B_Ik.dot(qid)
    return q

def q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik):
    q = th.mm(A_Ik,p) + th.mm(B_Ik,qid)
    return q

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
        jacobian[i] = inputs.grad.data
    return th.transpose(jacobian, dim0=0, dim1=1)


# ------------------ MAIN ----------------------

file_address = compas_tno.get('test.json')
# file_address = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form = FormDiagram.from_json(file_address)
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

# Compute the q's as tensor

q_ = q_from_qid_matrix_pythorch(qid_th, p_th, A_Ik_th, B_Ik_th)

# Compute Jacobian by mark

jac = compute_jacobian(qid_th.flatten, q_.flatten)
