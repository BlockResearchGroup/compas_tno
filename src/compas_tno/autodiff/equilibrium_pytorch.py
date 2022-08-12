from torch import mm
from torch import diagflat
from torch import tensor
from torch import cat
from torch import mul
from torch import div
from torch import solve
from torch import norm
from torch import zeros as thzeros
from torch.autograd.gradcheck import zero_gradients
from torch import float64

from numpy import hstack
from numpy import zeros
from numpy import array
from numpy import newaxis
from numpy import vstack
from numpy.linalg import pinv

from scipy.sparse import diags
from scipy.sparse import vstack as svstack
from scipy.sparse import csr_matrix

from compas.numerical import normrow
from compas.numerical import connectivity_matrix
from compas.utilities import geometric_key

from compas_tno.algorithms import find_independents
from compas_tno.algorithms import check_horizontal_loads
from compas_tno.algorithms import check_independents

import time


__all__ = [
    'q_from_variables_pytorch',
    'z_from_variables_pytorch',
    'reac_bound_variables_pytorch',
    'f_min_thrust_pytorch',
    'f_max_thrust_pytorch',
    'f_constraints_pytorch',
    'f_constraints_pytorch_MMA',
    'bounds_constraints_pytorch',
    'f_objective_pytorch',
    'compute_autograd',
    'compute_autograd_jacobian'
]


def q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep):
    k = len(ind)
    m = k + len(dep)
    q = thzeros(m, 1, dtype=float64)
    q[ind] = variables[:k]
    q[dep] = - Edinv_p_th + mm(EdinvEi_th, q[ind])
    return q


def z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci, Cit, Cf, pzfree, free, fixed):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    zfixed = variables[-Cf.shape[1]:]
    Q = diagflat(q)
    Ai = mm(mm(Cit, Q), Ci)
    Af = mm(mm(Cit, Q), Cf)
    b = pzfree - mm(Af, zfixed)
    zi, LU = solve(b, Ai)
    z = tensor(zeros((len(free)+len(fixed), 1)))
    z[free] = zi
    z[fixed] = zfixed
    return z


def reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, pfixed, xyz):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    zfixed = variables[-Cf.shape[1]:]
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    R = mm(CfQC, xyz) - pfixed
    length_x = mul(zfixed, abs(div(R[:, 0], R[:, 2])).reshape(-1, 1))  # Missing - s[0] in case the "datum is not at z=0"
    length_y = mul(zfixed, abs(div(R[:, 1], R[:, 2])).reshape(-1, 1))
    length = cat([length_x, length_y])
    return length


def f_min_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, xy, pfixed):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    Rh = mm(CfQC, xy) - pfixed[:, :2]
    R = norm(Rh, dim=1)
    f = sum(R)
    return f


def f_max_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, xy, pfixed):
    f = -1 * f_min_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, xy, pfixed)
    return f


def f_loadpath_pytorch(variables):

    return


def f_target_pytorch(variables):

    return


def f_objective_pytorch(variables, *args):
    Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, objective = args
    if objective == 'min':
        f = f_min_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, xy, pfixed)
    if objective == 'max':
        f = f_max_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, xy, pfixed)
    if objective == 'loadpath':
        raise NotImplementedError
    if objective == 'target':
        raise NotImplementedError
    return f


def f_constraints_pytorch(variables, *args):
    (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx,
     max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th) = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q[dep]
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        xyz = cat((xyz[:, :2], z), 1)
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = C_th[:, rol_x].t()
        Cfty = C_th[:, rol_y].t()
        rx_check = tensor(max_rol_rx) - abs(mm(Cftx, mm(U_th, q)) - tensor(px[rol_x]))
        ry_check = tensor(max_rol_ry) - abs(mm(Cfty, mm(V_th, q)) - tensor(py[rol_y]))
        constraints = cat([constraints, rx_check, ry_check])
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        A_q = mm(tensor(Asym), variables)
        constraints = cat([constraints, A_q])
    return constraints


def f_constraints_pytorch_MMA(variables, *args):  # THe equality constraints are transformed in 2 inequalities (symmetry constraint)
    (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx,
     max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th) = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = Cf_th[:, rol_x].transpose()
        Cfty = Cf_th[:, rol_y].transpose()
        rx_check = max_rol_rx - abs(mm(Cftx, mm(U_th, q)) - px[rol_x])
        ry_check = max_rol_ry - abs(mm(Cfty, mm(V_th, q)) - py[rol_y])
        constraints = cat([constraints, rx_check, ry_check])
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        A_q = mm(tensor(Asym), variables)
        constraints = cat([constraints, A_q, -1 * A_q])
    return constraints


def bounds_constraints_pytorch(variables, *args):
    (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx, max_rol_ry, rol_x,
     rol_y, px, py, Asym, U_th, V_th) = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = C_th[:, rol_x].t()
        Cfty = C_th[:, rol_y].t()
        rx_check = tensor(max_rol_rx) - abs(mm(Cftx, mm(U_th, q)) - tensor(px[rol_x]))
        ry_check = tensor(max_rol_ry) - abs(mm(Cfty, mm(V_th, q)) - tensor(py[rol_y]))
        constraints = cat([constraints, rx_check, ry_check])
    cu = [10e10]*len(constraints)
    cl = [0.0]*len(constraints)
    # End of inequality constraints
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        A_q = mm(tensor(Asym), variables)
        cu = hstack([cu, [0.0]*len(A_q)])
        cl = hstack([cl, [0.0]*len(A_q)])
    return cu, cl


def compute_autograd(variables, f):
    f.backward(retain_graph=True)
    grad = variables.grad.data
    return grad


def compute_autograd_jacobian(inputs, outputs):
    d_otp = outputs.size()[0]
    d_inp = inputs.size()[0]
    jacobian = tensor(zeros((d_otp, d_inp)))
    grad_output = tensor(zeros((d_otp, 1)))
    for i in range(d_otp):
        zero_gradients(inputs)
        grad_output.zero_()
        grad_output[i] = 1
        outputs.backward(grad_output, retain_graph=True)
        jacobian[i] = inputs.grad.data.flatten()
    return jacobian


def initialise_problem_torch(form, indset=None, printout=False, find_inds=True, tol=0.001):
    """ Initialise the problem for a given Form Diagram and return the set of matrices and vectors to optimise.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    printout : bool, optional
        Control output. By default False.
    find_inds : bool, optional
        If True will calculate the independents. Default is True
    indset : list, optional
        Independent set to use. If empty the independents are calculated normally. By default None
    tol: float, optional
        Tolerance for check the independent edges equilibrium. By default 0.001

    Returns
    -------
    args
        List of matrices and vectors used to perform the optimisation.

    """

    import torch as th

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()
    uv_i = form.uv_index()

    # Vertices and edges

    n = form.number_of_vertices()
    m = len(list(form.edges_where({'_is_edge': True})))
    fixed = [k_i[key] for key in form.fixed()]
    rol = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    sym = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free = list(set(range(n)) - set(fixed) - set(rol))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for i in range(n):
        key = i_k[i]
        if form.vertex_attribute(key, 'lb', None):
            lb_ind.append(i)
            lb.append(form.vertex_attribute(key, 'lb'))
        if form.vertex_attribute(key, 'ub', None):
            ub_ind.append(i)
            ub.append(form.vertex_attribute(key, 'ub'))

    lb = array(lb)
    ub = array(ub)
    lb.shape = (len(lb), 1)
    ub.shape = (len(ub), 1)

    # Co-ordinates and loads

    xyz = th.tensor(zeros((n, 3)))
    x = th.tensor(zeros((n, 1)))
    y = th.tensor(zeros((n, 1)))
    z = th.tensor(zeros((n, 1)))
    s = th.tensor(zeros((n, 1)))
    px = th.tensor(zeros((n, 1)))
    py = th.tensor(zeros((n, 1)))
    pz = th.tensor(zeros((n, 1)))
    s = th.tensor(zeros((n, 1)))
    w = th.tensor(zeros((n, 1)))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i] = vertex.get('x')
        y[i] = vertex.get('y')
        z[i] = vertex.get('z')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i] = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0)  # weight used in case of fiting...

    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]

    # If Rols are set up

    rol_x = []
    rol_y = []

    for key in form.vertices_where({'rol_x': True}):
        rol_x.append(k_i[key])

    for key in form.vertices_where({'rol_y': True}):
        rol_y.append(k_i[key])

    free_x = list(set(free) - set(rol_x))
    free_y = list(set(free) - set(rol_y))

    # C and E matrices

    C = th.tensor(connectivity_matrix(edges, 'array'))
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cftx = C[:, rol_x].transpose
    Cfty = C[:, rol_y].transpose
    Ct = C.transpose
    Cit = Ci.transpose
    Citx = Ct[free_x, :]  # .toarray()
    City = Ct[free_y, :]  # .toarray()
    uvw = C.dot(xyz)
    U = diags(uvw[:, 0].flatten())
    V = diags(uvw[:, 1].flatten())
    E = svstack((Citx.dot(U), City.dot(V)), dtype='csr').toarray()
    print('Equilibrium Matrix Shape: ', E.shape)

    start_time = time.time()

    # Independent and dependent branches

    if find_inds:

        if indset:
            ind = []
            for u, v in form.edges_where({'_is_edge': True}):
                if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                    ind.append(uv_i[(u, v)])
        else:
            ind = find_independents(E)

        k = len(ind)
        dep = list(set(range(m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print('Shape Equilibrium Matrix: ', E.shape)
            print('Found {0} independents'.format(k))
            print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

        for u, v in form.edges_where({'_is_edge': True}):
            form.edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

        Edinv = -csr_matrix(pinv(E[:, dep]))
        Ei = E[:, ind]

    else:
        k = m
        ind = list(set(range(m)))
        dep = []
        Edinv = []
        Ei = []

    # Set-up

    lh = normrow(C.dot(xy))**2
    p = vstack([px[free_x], py[free_y]])
    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]

    if any(p):
        check_hor = check_horizontal_loads(E, p)
        if check_hor:
            print('Horizontal Loads can be taken!')
        else:
            print('Horizontal Loads are not suitable for this FD!')

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k,
            lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty)
    args_inds = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Citx, City, Cf, U, V, p, px, py, pz, z, free_x, free_y, fixed, lh, sym, k)
    if find_inds is True:
        checked = check_independents(args_inds, tol=tol)
        if checked:
            print('Independents checked!')
            pass
        else:
            print('Warning: independent edges not equilibrated')

    return args
