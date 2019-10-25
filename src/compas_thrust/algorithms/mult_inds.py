# from numpy import abs
# from numpy import argmin
# from numpy import array
# from numpy import float64
# from numpy import dot
# from numpy import hstack
# from numpy import isnan
# from numpy import max
# from numpy import min
# from numpy import newaxis
# from numpy import sqrt
# from numpy import sum
# from numpy import vstack
# from numpy import zeros
# from numpy import ones
# from numpy import append
# from numpy.linalg import pinv
# from numpy.linalg import matrix_rank
# from numpy.random import rand
# from numpy.random import randint

# from scipy.linalg import svd
# from scipy.optimize import fmin_slsqp
# from scipy.sparse import csr_matrix
# from scipy.sparse import diags
# from scipy.sparse.linalg import spsolve
# import time

# from compas_tna.diagrams import FormDiagram

# from compas.numerical import connectivity_matrix
# from compas.numerical import devo_numpy
# from compas.numerical import equilibrium_matrix
# from compas.numerical import normrow
# from compas.numerical import nonpivots
# from compas.utilities import geometric_key

# from compas.geometry import intersection_line_line
# from compas.geometry import is_point_on_segment

# from compas_thrust.plotters.plotters import plot_form

# from compas_thrust.algorithms.equilibrium import zlq_from_qid
# from compas_thrust.algorithms.equilibrium import reactions
# from compas_thrust.algorithms.equilibrium import xlq_from_qid
# from compas_thrust.algorithms.equilibrium import ylq_from_qid

# from multiprocessing import Pool
# from random import shuffle
# from copy import deepcopy


# __author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
# __copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
# __license__   = 'MIT License'
# __email__     = 'mricardo@ethz.ch'


# __all__ = [
#     'optimise_single',
#     'optimise_multi',
#     'run_replicate',
#     'find_independents',
#     'initialize_problem',
# ]


# def optimise_single(form, solver='devo', polish='slsqp', qmin=1e-6, qmax=10, population=300, generations=500,
#                     printout=10, tol=0.01, plot=False, frange=[], indset=None, tension=False, planar=False,
#                     opt_max=False, t=None, random_ind=False, target=False, use_bounds=None, bounds_width = 5.0,
#                     summary = True, objective='loadpath', buttress = False):

#     """ Finds the optimised load-path for a FormDiagram.

#     Parameters
#     ----------
#     form : obj
#         The FormDiagram.
#     solver : str
#         Differential Evolution 'devo' or Genetic Algorithm 'ga' evolutionary solver to use.
#     polish : str
#         'slsqp' polish or None.
#     qmin : float
#         Minimum qid value.
#     qmax : float
#         Maximum qid value.
#     population : int
#         Number of agents for the evolution solver.
#     generations : int
#         Number of generations for the evolution solver.
#     printout : int
#         Frequency of print output to the terminal.
#     tol : float
#         Tolerance on horizontal force balance.
#     plot : bool
#         Plot progress of the evolution.
#     frange : list
#         Minimum and maximum function value to plot.
#     indset : list
#         Independent set to use.
#     tension : bool
#         Allow tension edge force densities (experimental).
#     planar : bool
#         Only consider the x-y plane.

#     Returns
#     -------
#     float
#         Optimum load-path value.
#     list
#         Optimum qids

#     """

#     if printout:
#         print('\n' + '-' * 50)
#         print('Load-path optimisation started')
#         print('-' * 50)

#     # Mapping

#     k_i  = form.key_index()
#     i_k  = form.index_key()
#     i_uv = form.index_uv()
#     uv_i = form.uv_index()

#     # Vertices and edges

#     n     = form.number_of_vertices()
#     m     = form.number_of_edges()
#     fixed = [k_i[key] for key in form.fixed()]
#     rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
#     edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
#     sym   = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
#     free  = list(set(range(n)) - set(fixed) - set(rol))

#     # Constraints

#     lb_ind = []
#     ub_ind = []
#     lb = []
#     ub = []
#     for key, vertex in form.vertex.items():
#         if vertex.get('lb', None):
#             lb_ind.append(k_i[key])
#             lb.append(vertex['lb'])
#         if vertex.get('ub', None):
#             ub_ind.append(k_i[key])
#             ub.append(vertex['ub'])
#     lb = array(lb)
#     ub = array(ub)
#     lb.shape = (len(lb),1)
#     ub.shape = (len(ub),1)

#     # Co-ordinates and loads

#     xyz = zeros((n, 3))
#     x   = zeros((n, 1))
#     y   = zeros((n, 1))
#     z   = zeros((n, 1))
#     s   = zeros((n, 1))
#     px  = zeros((n, 1))
#     py  = zeros((n, 1))
#     pz  = zeros((n, 1))
#     w   = zeros((n, 1))

#     for key, vertex in form.vertex.items():
#         i = k_i[key]
#         xyz[i, :] = form.vertex_coordinates(key)
#         x[i]  = vertex.get('x')
#         y[i]  = vertex.get('y')
#         z[i] = vertex.get('z')
#         px[i] = vertex.get('px', 0)
#         py[i] = vertex.get('py', 0)
#         pz[i] = vertex.get('pz', 0)
#         s[i]  = vertex.get('target', 0)
#         w[i] = vertex.get('weight', 1.0)

#     Wfree = diags(w[free].flatten())
#     if printout:
#         print('Max Z is', max(z))
#     xy = xyz[:, [0,1]]
#     xz = xyz[:, [0,2]]
#     yz = xyz[:, [1,2]]
#     print(xyz[0])
#     print(xy[0])
#     print(xz[0])
#     print(yz[0])

#     # Temporary anchors location

#     anchors = []

#     # C and E matrices

#     C   = connectivity_matrix(edges, 'csr')
#     Ci  = C[:, free]
#     Cf  = C[:, fixed]
#     Ct = C.transpose()
#     Cit = Ci.transpose()
#     Exy   = equilibrium_matrix_3D(C, xyz, free, 'csr', 'xy').toarray()
#     Exz   = equilibrium_matrix_3D(C, xyz, free, 'csr', 'xz').toarray()
#     Eyz   = equilibrium_matrix_3D(C, xyz, free, 'csr', 'yz').toarray()
#     uvw = C.dot(xyz)
#     U   = uvw[:, 0]
#     V   = uvw[:, 1]
#     W   = uvw[:, 2]

#     # Independent and dependent branches

#     indz = find_independents(Exy)
#     indy = find_independents(Exz)
#     indx = find_independents(Eyz)

#     print('Shape Z Equilibrium Matrix: ', Exy.shape)
#     print('Rank Equilibrium Matrix: ', matrix_rank(Exy))
#     print('Indeoendents Z: {0} '.format(len(indz)))
#     print('Shape Y Equilibrium Matrix: ', Exz.shape)
#     print('Rank Equilibrium Matrix: ', matrix_rank(Exz))
#     print('Indeoendents Y: {0}'.format(len(indy)))
#     print('Shape X Equilibrium Matrix: ', Eyz.shape)
#     print('Rank Equilibrium Matrix: ', matrix_rank(Eyz))
#     print('Indeoendents X: {0}'.format(len(indx)))

#     print(indz)
#     print(indy)
#     print(indx)
#     # ind = indz
#     kz   = len(indz)
#     ky   = len(indy)
#     kx   = len(indx)
#     depz = list(set(range(m)) - set(indz))
#     depy = list(set(range(m)) - set(indy))
#     depx = list(set(range(m)) - set(indx))

#     for u, v in form.edges():
#         form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in indz else False)

#     # plot_form(form, fix_width=True, max_width = 2.0, show_q= False).show()

#     # for u, v in form.edges():
#     #     form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in indy else False)

#     # plot_form(form, fix_width=True, max_width = 2.0, show_q= False).show()

#     # for u, v in form.edges():
#     #     form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in indx else False)

#     # plot_form(form, fix_width=True, max_width = 2.0, show_q= False).show()

#     if printout:
#         print('Shell element has {0}, {1}, {2} (Z, Y, X) independent branches '.format(len(indz),len(indy),len(indx)))
#         print('Shape Equilibrium Matrices: ', Exy.shape, Exz.shape, Eyz.shape)
#         print('Rank Equilibrium Matrix: ', matrix_rank(Exy), matrix_rank(Exz), matrix_rank(Eyz))

#     # Set-up

#     lxy     = normrow(C.dot(xy))**2
#     lxz     = normrow(C.dot(xz))**2
#     lyz     = normrow(C.dot(yz))**2
#     Edinv_z  = -csr_matrix(pinv(Exy[:, depz]))
#     Edinv_y  = -csr_matrix(pinv(Exz[:, depy]))
#     Edinv_x  = -csr_matrix(pinv(Eyz[:, depx]))
#     Eiz     = Exy[:, indz]
#     Eiy     = Exz[:, indy]
#     Eix     = Eyz[:, indx]
#     p_z     = vstack([px[free], py[free]])
#     p_y     = vstack([px[free], pz[free]])
#     p_x     = vstack([py[free], pz[free]])
#     q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

#     print('Shape Inverse Matrices: ', Edinv_z.shape, Edinv_y.shape, Edinv_x.shape)
#     print('Rank Inverses Matrix: ', matrix_rank(Edinv_z), matrix_rank(Edinv_y), matrix_rank(Edinv_x))

#     if objective == 'loadpath':
#         t = None

#     if buttress:
#         b = {}
#         for key in form.vertices_where({'is_fixed': True}):
#             try:
#                 b[k_i[key]] = form.get_vertex_attributes(key, 'b')
#             except:
#                 pass
#                 print(b)
#     else:
#         b = None

#     if printout and buttress:
#         print('Butress length {0}'.format(b))
    
#     try:
#         joints = form.attributes['joints']
#     except:
#         joints = None

#     if printout and joints:
#         print('Joins Data', joints)

#     args = (q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x)

#     if use_bounds:
#         bounds = []
#         for i in ind:
#             u,v = i_uv[i]
#             qi = form.get_edge_attribute((u,v), 'q')
#             if form.get_edge_attribute((u,v), 'is_symmetry') is True:
#                 bounds.append([qmin, qmax])
#             else:
#                 if qi > bounds_width + qmin:
#                     bounds.append([qi-bounds_width,qi+bounds_width])
#                 else:
#                     bounds.append([qmin,qi+bounds_width])
#     else:
#         boundsz = [[qmin, qmax]] * kz
#         boundsy = [[qmin, qmax]] * ky
#         boundsx = [[qmin, qmax]] * kx
#     if t:
#         bounds.append([-t,t])
    

#     # Horizontal checks

#     checked = True

#     _q = deepcopy(q)

#     if tol > 0:
#         for i in range(10**3):
#             _q[indz, 0] = rand(kz) * qmax
#             _q[depz] = -Edinv_z.dot(p_z - Eiz.dot(_q[indz]))
#             Rx = Cit.dot(U * _q.ravel()) - px[free].ravel()
#             Ry = Cit.dot(V * _q.ravel()) - py[free].ravel()
#             Rz = Cit.dot(W * _q.ravel())
#             R  = max(sqrt(Rx**2 + Ry**2))
#             # print(sum(Rz))
#             if R > tol:
#                 checked = False
#                 print('Failed checks on Z. Iteration {0}'.format(i))
#                 print('Rx: {0} Rz: {1}'.format(sum(abs(Rx)),sum(abs(Ry))))
#                 break

#     _q = deepcopy(q)

#     if tol > 0:
#         for i in range(10**3):
#             _q[indy, 0] = rand(ky) * qmax
#             _q[depy] = -Edinv_y.dot(p_y - Eiy.dot(_q[indy]))
#             Rx = Cit.dot(U * _q.ravel()) - px[free].ravel()
#             Rz = Cit.dot(W * _q.ravel()) - pz[free].ravel()
#             R  = max(sqrt(Rx**2 + Rz**2))
#             # print(sum(Rz))
#             if R > tol:
#                 checked = False
#                 # print('\nFailed checks on Y. Iteration {0}'.format(i))
#                 # print('Rx: {0} Rz: {1}'.format(sum(abs(Rx)),sum(abs(Rz))))
#                 # break

#     _q = deepcopy(q)

#     if tol > 0:
#         for i in range(5):
#             _q[indx, 0] = rand(kx) * qmax
#             _q[depx] = -Edinv_x.dot(p_x - Eix.dot(_q[indx]))
#             Rz = Cit.dot(W * _q.ravel()) - pz[free].ravel()
#             Ry = Cit.dot(V * _q.ravel()) - py[free].ravel()
#             R  = max(sqrt(Rz**2 + Ry**2))
#             # print(sum(Rz))
#             if R > tol:
#                 checked = False
#                 # print('\nFailed checks on X. Iteration {0}'.format(i))
#                 # print('Rz: {0} Ry: {1}\n'.format(sum(abs(Rz)),sum(abs(Ry))))
#                 # break

#     if checked:

#         plot_form(form).show()

#         # Define Objective

#         if objective=='loadpath':
#             fdevo, fslsqp, fieq = _fint, _fint_, _fieq
#         if objective=='target':
#             fdevo, fslsqp, fieq = fbf, fbf_dome_, _fieq
#         if objective == 'min':
#             fdevo, fslsqp, fieq = fmin, fmin, _fieq_bounds
#         if objective=='max':
#             fdevo, fslsqp, fieq = fmax, fmax, _fieq_bounds
#         if objective=='bounds':
#             fdevo, fslsqp, fieq = fbounds, fbounds, _fieq_bounds

#         # Solve with gradient in Z, X, Y .
        
#         # Solve in Z:

#         qi = q[indz]
#         if t:
#             qopt = append(qi,t)
#         else:
#             qopt = qi.reshape(kz,)
#         fopt = fslsqp(qopt,*args)
#         if printout:
#             print('Gradient Method: Initial Value: {0}'.format(fopt))

#         if polish == 'slsqp':
#             fopt_, qopt_ = _slsqp(fslsqp, qopt, boundsz, printout, fieq, args)
#             if t:
#                 qopt_, t_slsqp = qopt_[:kz], qopt_[-1]
#                 qopt, t = qopt[:kz], qopt[-1]
#             q1 = zlq_from_qid(qopt_, args)[2]
#             if fopt_ < fopt:
#                 if (min(q1) > -0.001 and not tension) or tension:
#                     fopt, qopt = fopt_, qopt_
#                     if t:
#                         t = t_slsqp
#         if t:
#             form.attributes['offset'] = t

#         z, _, q, q_ = zlq_from_qid(qopt, args)

#         for i in range(n):
#             key = i_k[i]
#             form.set_vertex_attribute(key=key, name='z', value=float(z[i]))
#             form.set_vertex_attribute(key=key, name='y', value=float(y[i]))
#             form.set_vertex_attribute(key=key, name='x', value=float(x[i]))

#         for c, qi in enumerate(list(q_.ravel())):
#             u, v = i_uv[c]
#             form.set_edge_attribute((u, v), 'q', float(qi))

#         plot_form(form).show()

#         # Set up for Y

#         for key, vertex in form.vertex.items():
#             i = k_i[key]
#             xyz[i, :] = form.vertex_coordinates(key)
#         xz = xyz[:, [0,2]]
#         Exz   = equilibrium_matrix_3D(C, xyz, free, 'csr', 'xz').toarray()
#         lxz     = normrow(C.dot(xz))**2
#         Edinv_y  = -csr_matrix(pinv(Exz[:, depy]))
#         Eiy     = Exz[:, indy]
#         args =  (q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x)
 
#         # Solve in Y:

#         qi = q[indy]
#         if t:
#             qopt = append(qi,t)
#         else:
#             qopt = qi.reshape(ky,)
#         fopt = fbf_dome_y(qopt,*args)
#         if printout:
#             print('Gradient Method for Y: Initial Value: {0}'.format(fopt))

#         if polish == 'slsqp':
#             fopt_, qopt_ = _slsqp(fbf_dome_y, qopt, boundsy, printout, _fieqy, args)
#             if t:
#                 qopt_, t_slsqp = qopt_[:ky], qopt_[-1]
#                 qopt, t = qopt[:ky], qopt[-1]
#             q1 = ylq_from_qid(qopt_, args)[2]
#             if fopt_ < fopt:
#                 if (min(q1) > -0.001 and not tension) or tension:
#                     fopt, qopt = fopt_, qopt_
#                     if t:
#                         t = t_slsqp
#         if t:
#             form.attributes['offset'] = t

#         y, _, q, q_ = ylq_from_qid(qopt, args)

#         for i in range(n):
#             key = i_k[i]
#             form.set_vertex_attribute(key=key, name='z', value=float(z[i]))
#             form.set_vertex_attribute(key=key, name='y', value=float(y[i]))
#             form.set_vertex_attribute(key=key, name='x', value=float(x[i]))

#         for c, qi in enumerate(list(q_.ravel())):
#             u, v = i_uv[c]
#             form.set_edge_attribute((u, v), 'q', float(qi))

#         plot_form(form).show()

#         # Set up for X

#         for key, vertex in form.vertex.items():
#             i = k_i[key]
#             xyz[i, :] = form.vertex_coordinates(key)
#         yz = xyz[:, [1,2]]
#         Eyz   = equilibrium_matrix_3D(C, xyz, free, 'csr', 'yz').toarray()
#         lyz     = normrow(C.dot(yz))**2
#         Edinv_x  = -csr_matrix(pinv(Eyz[:, depx]))
#         Eix     = Eyz[:, indx]
#         args =  (q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x)
#         # Solve in X:

#         qi = q[indx]
#         if t:
#             qopt = append(qi,t)
#         else:
#             qopt = qi.reshape(kx,)
#         fopt = fbf_dome_x(qopt,*args)
#         if printout:
#             print('Gradient Method for X: Initial Value: {0}'.format(fopt))

#         if polish == 'slsqp':
#             fopt_, qopt_ = _slsqp(fbf_dome_x, qopt, boundsx, printout, _fieqx, args)
#             if t:
#                 qopt_, t_slsqp = qopt_[:kx], qopt_[-1]
#                 qopt, t = qopt[:kx], qopt[-1]
#             q1 = xlq_from_qid(qopt_, args)[2]
#             if fopt_ < fopt:
#                 if (min(q1) > -0.001 and not tension) or tension:
#                     fopt, qopt = fopt_, qopt_
#                     if t:
#                         t = t_slsqp
#         if t:
#             form.attributes['offset'] = t

#         for i in range(n):
#             key = i_k[i]
#             form.set_vertex_attribute(key=key, name='z', value=float(z[i]))
#             form.set_vertex_attribute(key=key, name='y', value=float(y[i]))
#             form.set_vertex_attribute(key=key, name='x', value=float(x[i]))

#         for c, qi in enumerate(list(q_.ravel())):
#             u, v = i_uv[c]
#             form.set_edge_attribute((u, v), 'q', float(qi))

#         x, _, q, q_ = xlq_from_qid(qopt, args)

#         plot_form(form).show()

#         # Update FormDiagram

#         for i in range(n):
#             key = i_k[i]
#             form.set_vertex_attribute(key=key, name='z', value=float(z[i]))
#             form.set_vertex_attribute(key=key, name='y', value=float(y[i]))
#             form.set_vertex_attribute(key=key, name='x', value=float(x[i]))

#         plot_form(form).show()

#         for c, qi in enumerate(list(q_.ravel())):
#             u, v = i_uv[c]
#             form.set_edge_attribute((u, v), 'q', float(qi))

#         # Relax

#         # q    = array([attr['q'] for u, v, attr in form.edges(True)])
#         # Q    = diags(q)
#         # CitQ = Cit.dot(Q)
#         # Di   = CitQ.dot(Ci)
#         # Df   = CitQ.dot(Cf)
#         # bx   = px - Df.dot(x[fixed])
#         # by   = py - Df.dot(y[fixed])
#         # # bz   = pz - Df.dot(z[fixed])
#         # x[free, 0] = spsolve(Di, bx)
#         # y[free, 0] = spsolve(Di, by)
#         # # z[free, 0] = spsolve(Di, bz)

#         # for i in range(n):
#         #     form.set_vertex_attributes(key=i_k[i], names='xyz', values=[float(j) for j in [x[i], y[i], z[i]]])
        
#         lp = 0
#         for u, v in form.edges():
#             if form.get_edge_attribute((u, v), 'is_symmetry') is False:
#                 qi = form.get_edge_attribute((u, v), 'q')
#                 li = form.edge_length(u, v)
#                 lp += abs(qi) * li**2

#         form.attributes['loadpath'] = lp

#         if objective == 'loadpath':
#             fopt = lp

#         CfQ = Cf.transpose().dot(diags(q.flatten()))

#         if objective == 'min':
#             fopt = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed]).flatten()[0]
#             [[fopt]] = fopt

#         if objective == 'max':
#             fopt = 1 / (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])[0][0]
#             [[fopt]] = fopt

#         reactions(form, plot = plot)

#         if summary:
#             print('\n' + '-' * 50)
#             print('qid range : {0:.3f} : {1:.3f}'.format(min(qopt), max(qopt)))
#             print('q range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))
#             print('fopt      : {0:.3f}'.format(fopt))
#             if t:
#                 print('t      : {0:.3f}'.format(t))
#             print('-' * 50 + '\n')

#         return fopt, qopt

#     else:

#         if printout:
#             print('Horizontal equillibrium checks failed')

#         return None, None


# def equilibrium_matrix_3D(C, xyz, free, rtype='array', axis='xy'):

#     from scipy.sparse import vstack as svstack
#     from numpy import asarray

#     xyz = asarray(xyz, dtype=float)
#     C   = csr_matrix(C)
#     uvw  = C.dot(xyz)
#     U   = diags([uvw[:, 0].flatten()], [0])
#     V   = diags([uvw[:, 1].flatten()], [0])
#     W   = diags([uvw[:, 2].flatten()], [0])
#     Ct  = C.transpose()
#     Cti = Ct[free, :]
#     if axis == 'xy':
#         E   = svstack((Cti.dot(U), Cti.dot(V)))
#     if axis == 'xz':
#         E   = svstack((Cti.dot(U), Cti.dot(W)))
#     if axis == 'yz':
#         E   = svstack((Cti.dot(V), Cti.dot(W)))
#     return _return_matrix(E, rtype)

# def _return_matrix(M, rtype):
#     if rtype == 'list':
#         return M.toarray().tolist()
#     if rtype == 'array':
#         return M.toarray()
#     if rtype == 'csr':
#         return M.tocsr()
#     if rtype == 'csc':
#         return M.tocsc()
#     if rtype == 'coo':
#         return M.tocoo()
#     return M

# def _fint(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     f = dot(abs(q.transpose()), l2)

#     if isnan(f):
#         return 10**10

#     else:

#         if not tension:
#             f += sum((q[q < 0] - 5)**4)

#         Rx = Cit.dot(U * q_.ravel()) - px[free].ravel()
#         Ry = Cit.dot(V * q_.ravel()) - py[free].ravel()
#         Rh = Rx**2 + Ry**2
#         Rm = max(sqrt(Rh))
#         if Rm > tol:
#             f += sum(Rh - tol + 5)**4

#         return f

# def fmin(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     CfQ = Cf.transpose().dot(diags(q_.flatten()))
#     f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])
#     # f +=  pz.transpose().dot(z[free])/1000
#     # print( (CfQ.dot(U[:,newaxis])) )
#     # print( (CfQ.dot(V[:,newaxis])) )
#     # f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

#     if isnan(f) == True or any(qid) == False:
#         return 10**10
#     else:

#         if not tension:
#             f += sum((q[q < 0] - 10)**4)
#         Rx = Ct.dot(U * q_.ravel()) - px.ravel()
#         Ry = Ct.dot(V * q_.ravel()) - py.ravel()
#         Rh = Rx**2 + Ry**2
#         Rm = max(sqrt(Rh[free]))
        
#         if Rm > tol:
#             f += sum(Rh - tol + 5)**4

#         if lb_ind:
#             if t is not 0.0:
#                 lb = lb - t
#             z_lb    = z[lb_ind]
#             log_lb  = z_lb < lb
#             diff_lb = z_lb[log_lb] - lb[log_lb]
#             pen_lb  = sum(abs(diff_lb) + 5)**4
#             f += pen_lb

#         if ub_ind:
#             if t is not 0.0:
#                 ub = ub - t
#             z_ub    = z[ub_ind]
#             log_ub  = z_ub > ub
#             diff_ub = z_ub[log_ub] - ub[log_ub]
#             pen_ub  = sum(abs(diff_ub) + 5)**4
#             f += pen_ub

#         if b:
#             W = C.dot(z)[:,0]
#             Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
#             for key in b:
#                 scl = t/Rz[key]
#                 x_comp = abs(scl * Rx[key])
#                 y_comp = abs(scl * Ry[key])
#                 if x_comp > abs(b[key][0][0]):
#                     f += (x_comp - b[key][0][0] + 5) **2
#                 if y_comp > abs(b[key][0][1]):
#                     f += (y_comp - b[key][0][1] + 5) **2

#         if joints:
#             for jt in joints:
#                 limit = [jt[0],jt[1]]
#                 u, v = k_i[jt[2][0]] , k_i[jt[2][1]]
#                 x_ = list(x)
#                 y_ = list(y)
#                 z_ = list(z)
#                 thrust = [[x_[u],y_[u],z_[u]+t],[x_[v],y_[v],z_[v]+t]]
#                 pt = intersection_line_line(limit,thrust)[0]
#                 if pt == None or is_point_on_segment(pt,limit,tol=1e-6) == False:
#                     f += 10 ** 5
                
#         return f

# def _fmin(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args

#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     CfQ = Cf.transpose().dot(diags(q.flatten()))
#     f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

#     if isnan(f) == True or any(qid) == False:
#         return 10**10

#     return f

# def fbounds(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args
#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     f = 0

#     if isnan(f):
#         return 10**10

#     else:

#         if not tension:
#             f += sum((q[q < 0] - 5)**4)

#         Rx = Cit.dot(U * q_.ravel()) - px[free].ravel()
#         Ry = Cit.dot(V * q_.ravel()) - py[free].ravel()
#         Rh = Rx**2 + Ry**2
#         Rm = max(sqrt(Rh))
#         if Rm > tol:
#             f += sum(Rh - tol + 5)**4

#         if lb_ind:
#             if t is not 0.0:
#                 lb = lb - t
#             z_lb    = z[lb_ind]
#             log_lb  = z_lb < lb
#             diff_lb = z_lb[log_lb] - lb[log_lb]
#             pen_lb  = sum(abs(diff_lb) + 5)**4
#             f += pen_lb

#         if ub_ind:
#             if t is not 0.0:
#                 ub = ub - t
#             z_ub    = z[ub_ind]
#             log_ub  = z_ub > ub
#             diff_ub = z_ub[log_ub] - ub[log_ub]
#             pen_ub  = sum(abs(diff_ub) + 5)**4
#             f += pen_ub

#         if b:
#             for key in b:
#                 scl = t/Rz[key]
#                 x_comp = abs(scl * Rx[key])
#                 y_comp = abs(scl * Ry[key])
#                 # print(b[key],scl,x_comp, y_comp)
#                 if x_comp > abs(b[key][0][0]):
#                     f += (x_comp - b[key][0][0] + 5) **2
#                 if y_comp > abs(b[key][0][1]):
#                     f += (y_comp - b[key][0][1] + 5) **2

#         return f

# def fmax(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args
    
#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     CfQ = Cf.transpose().dot(diags(q_.flatten()))
#     f = 1 / ((CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed]))

#     if isnan(f) == True or any(qid) == False:
#         return 10**10

#     else:

#         if not tension:
#             f += sum((q[q < 0] - 5)**4)
#         W = C.dot(z)[:,0]
#         Rx = Ct.dot(U * q_.ravel()) - px.ravel()
#         Ry = Ct.dot(V * q_.ravel()) - py.ravel()
#         Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
#         Rh = Rx**2 + Ry**2
#         Rm = max(sqrt(Rh[free]))
#         if Rm > tol:
#             f += sum(Rh - tol + 5)**4

#         if lb_ind:
#             if t is not 0.0:
#                 lb = lb - t
#             z_lb    = z[lb_ind]
#             log_lb  = z_lb < lb
#             diff_lb = z_lb[log_lb] - lb[log_lb]
#             pen_lb  = sum(abs(diff_lb) + 2)**4
#             f += pen_lb

#         if ub_ind:
#             if t is not 0.0:
#                 ub = ub - t
#             z_ub    = z[ub_ind]
#             log_ub  = z_ub > ub
#             diff_ub = z_ub[log_ub] - ub[log_ub]
#             pen_ub  = sum(abs(diff_ub) + 2)**4
#             f += pen_ub

#         if b:
#             for key in b:
#                 scl = t/Rz[key]
#                 x_comp = abs(scl * Rx[key])
#                 y_comp = abs(scl * Ry[key])
#                 # print(b[key],scl,x_comp, y_comp)
#                 if x_comp > abs(b[key][0][0]):
#                     f += (x_comp - b[key][0][0] + 5) **2
#                 if y_comp > abs(b[key][0][1]):
#                     f += (y_comp - b[key][0][1] + 5) **2

#         return f

# def fbf(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     z = z + t
#     z_s = abs(z[free] - s[free])
#     f = Wfree.dot(z_s.transpose().dot(z_s))
#     f = f[0,0]

#     if not tension:
#             f += sum((q[q < 0] - 5)**4)

#     if isnan(f):
#         return 10**10

#     else:
#         return f

# def _fint_(qid, *args):

#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     f = dot(abs(q.transpose()), l2)

#     if isnan(f):
#         return 10**10

#     return f

# def fbf_(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

#     qid, t = qid[:k], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     z_s = abs(z[free] - s[free])
#     f = z_s.transpose().dot(z_s)
#     f = f[0,0]

#     if isnan(f):
#         return 10**10

#     else:
#         return f

# def fbf_dome_(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args

#     qid, t = qid[:kz], qid[-1]
#     z, l2, q, q_ = zlq_from_qid(qid, args)
#     exp = sqrt(x[free]**2 + y[free]**2 + z[free]**2)
#     f = sum(abs(exp - 10.0 * ones((len(free),)))) # 10 is the radius

#     if isnan(f):
#         return 10**10

#     else:
#         return f

# def fbf_dome_y(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args
    
#     qid, t = qid[:len(indy)], qid[-1]
#     y, l2, q, q_ = ylq_from_qid(qid, args)
#     exp = sqrt(x[free]**2 + y[free]**2 + z[free]**2)
#     f = sum(abs(exp - 10.0 * ones((len(free),)))) # 10 is the radius

#     if isnan(f):
#         return 10**10

#     else:
#         return f

# def fbf_dome_x(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args
#     qid, t = qid[:len(indx)], qid[-1]
#     x, l2, q, q_ = xlq_from_qid(qid, args)
#     exp = sqrt(x[free]**2 + y[free]**2 + z[free]**2)
#     f = sum(abs(exp - 10.0 * ones((len(free),)))) # 10 is the radius

#     if isnan(f):
#         return 10**10

#     else:
#         return f

# def _fieq(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args
    
#     q[indz, 0], t = qid[:kz], qid[-1]
#     q[depz] = -Edinv_z.dot(p_z - Eiz.dot(q[indz]))

#     if not tension:
#         return hstack([q.ravel() + 10**(-5)])
#         # return hstack([Rm + pen_lb + pen_ub])
#     return [tol - Rm - pen_lb - pen_ub]

# def _fieqy(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args
    
#     q[indy, 0], t = qid[:len(indy)], qid[-1]
#     # z, l2, q, q_ = zlq_from_qid(qid, args) # Added
#     depy = list(set(range(len(q))) - set(indy))
#     q[depy] = -Edinv_y.dot(p_y - Eiy.dot(q[indy]))

#     if not tension:
#         return hstack([q.ravel() + 10**(-5)])
#         # return hstack([Rm + pen_lb + pen_ub])
#     return [tol - Rm - pen_lb - pen_ub]

# def _fieqx(qid, *args):

#     q, indz, depz, Edinv_z, Eiz, C, Ct, Ci, Cit, Cf, U, V, p_z, px, py, pz, tol, z, free, fixed, planar, lxy, sym, tension, kz, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i, indy, Edinv_y, Eiy, p_y, lxz, y, indx, Edinv_x, Eix, p_x, lyz, W, x = args
    
#     q[indx, 0], t = qid[:len(indx)], qid[-1]
#     # z, l2, q, q_ = zlq_from_qid(qid, args) # Added
#     depx = list(set(range(len(q))) - set(indx))
#     q[depx] = -Edinv_x.dot(p_x - Eiy.dot(q[indx]))

#     if not tension:
#         return hstack([q.ravel() + 10**(-5)])
#         # return hstack([Rm + pen_lb + pen_ub])
#     return [tol - Rm - pen_lb - pen_ub]

# def _fieq_ub_lb(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

#     q[ind, 0], t = qid[:k], qid[-1]
#     q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
#     q_ = 1 * q
#     q[sym] *= 0

#     Rm = 0

#     pen_lb = 0.0
#     pen_ub = 0.0

#     if lb_ind:
#         if t is not 0.0:
#             lb = lb - t
#         z_lb    = z[lb_ind]
#         log_lb  = z_lb < lb
#         diff_lb = z_lb[log_lb] - lb[log_lb]
#         pen_lb  = sum(abs(diff_lb) + 2)**4

#     if ub_ind:
#         if t is not 0.0:
#             ub = ub - t
#         z_ub    = z[ub_ind]
#         log_ub  = z_ub > ub
#         diff_ub = z_ub[log_ub] - ub[log_ub]
#         pen_ub  = sum(abs(diff_ub) + 2)**4

#     if not tension:
#         return hstack([pen_ub + pen_ub])
#         # return hstack([Rm + pen_lb + pen_ub])
#     return [tol - Rm - pen_lb - pen_ub]

# def _fieq_bounds(qid, *args):

#     q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args
    
#     q[ind, 0], t = qid[:k], qid[-1]
#     q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
#     q_ = 1 * q
#     q[sym] *= 0
    
#     rbut = 0.0

#     W = C.dot(z)[:,0]
#     Rx = Ct.dot(U * q_.ravel()) - px.ravel()
#     Ry = Ct.dot(V * q_.ravel()) - py.ravel()
#     Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
#     Rh = Rx**2 + Ry**2
#     Rm = max(sqrt(Rh[free]))

#     pen_lb = 0.0
#     pen_ub = 0.0
#     joint = 0.0

#     if lb_ind:
#         if t is not 0.0:
#             lb = lb - t
#         z_lb    = z[lb_ind]
#         log_lb  = z_lb < lb
#         diff_lb = z_lb[log_lb] - lb[log_lb]
#         pen_lb  = sum(abs(diff_lb) + 2)**4

#     if ub_ind:
#         if t is not 0.0:
#             ub = ub - t
#         z_ub    = z[ub_ind]
#         log_ub  = z_ub > ub
#         diff_ub = z_ub[log_ub] - ub[log_ub]
#         pen_ub  = sum(abs(diff_ub) + 2)**4

#     if b:
#             for key in b:
#                 scl = t/Rz[key]
#                 x_comp = abs(scl * Rx[key])
#                 y_comp = abs(scl * Ry[key])
#                 if x_comp > abs(b[key][0][0]):
#                     rbut += (x_comp - b[key][0][0] + 5) **2
#                 if y_comp > abs(b[key][0][1]):
#                     rbut += (y_comp - b[key][0][1] + 5) **2

#     if joints:
#         for jt in joints:
#             limit = [jt[0],jt[1]]
#             u, v = k_i[jt[2][0]] , k_i[jt[2][1]]
#             x_ = list(x)
#             y_ = list(y)
#             z_ = list(z)
#             thrust = [[x_[u],y_[u],z_[u]+t],[x_[v],y_[v],z_[v]+t]]
#             pt = intersection_line_line(limit,thrust)[0]
#             if pt == None or is_point_on_segment(pt,limit,tol=1e-6) == False:
#                 joint += 10 ** 5

#     if not tension:
#         return hstack([Rm + pen_lb + pen_ub])
#         # return hstack([q.ravel() + 10**(-5), tol - Rm - pen_lb - pen_ub])
#         # return hstack([Rm + pen_lb + pen_ub])
#     return [tol - Rm - pen_lb - pen_ub - rbut - joint]

# def _slsqp(fn, qid0, bounds, printout, fieq, args):

#     pout = 2 if printout else 0
#     opt  = fmin_slsqp(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=500, f_ieqcons=fieq)

#     return opt[1], opt[0]

# def _diff_evo(fn, bounds, population, generations, printout, plot, frange, args):

#     return devo_numpy(fn=fn, bounds=bounds, population=population, generations=generations, printout=printout,
#                       plot=plot, frange=frange, args=args)

# def _worker(data):

#     try:

#         i, form, save_figs, qmin, qmax, population, generations, simple, tension, tol = data
#         fopt, qopt = optimise_single(form, qmin=qmin, qmax=qmax, population=population, generations=generations,
#                                      printout=0, tension=tension, tol=tol)

#         if isnan(fopt):
#             fopt = 10**10

#         print('Trial: {0} - Optimum: {1:.1f}'.format(i, fopt))

#         if save_figs:
#             plotter = plot_form(form, radius=0.1, fix_width=False, simple=simple)
#             plotter.save('{0}trial_{1}-fopt_{2:.6f}.png'.format(save_figs, i, fopt))
#             del plotter

#         return fopt, form

#     except:

#         print('Trial: {0} - FAILED'.format(i))

#         return 10**10, None

# def optimise_multi(form, trials=10, save_figs='', qmin=0.001, qmax=5, population=300, generations=500, simple=False,
#                    tension=False, tol=0.001, opt_max=False):

#     """ Finds the optimised load-path for multiple d FormDiagrams.

#     Parameters
#     ----------
#     form : obj
#         FormDiagram to analyse.
#     trials : int
#         Number of trials to perform.
#     save_figs : str
#         Directory to save plots.
#     qmin : float
#         Minimum qid value.
#     qmax : float
#         Maximum qid value.
#     population : int
#         Number of agents for the evolution solver.
#     generations : int
#         Number of generations for the evolution solver.
#     simple : bool
#         Simple red and blue colour plotting.
#     tension : bool
#         Allow tension edge force densities (experimental).
#     tol : float
#         Tolerance on horizontal force balance.

#     Returns
#     -------
#     list
#         Optimum load-path for each trial.
#     list
#         Each resulting trial FormDiagram.
#     int
#         Index of the optimum.

#     """

#     data = [(i, _form(form), save_figs, qmin, qmax, population, generations, simple, tension, tol)
#             for i in range(trials)]

#     fopts, forms = zip(*Pool().map(_worker, data))
#     best = argmin(fopts)
#     if opt_max:
#         best = argmax(fopts)

#     print('Best: {0} - fopt {1:.1f}'.format(best, fopts[best]))

#     return fopts, forms, best

# def run_replicate(form, file_complete, delete_face=False, plots=False, save=True):

#     print('\nOptimization and Replicate of Form ')

#     form = FormDiagram.from_json(file)
#     form = _form(form)

#     plot_form(form, radius=0.05).show()

#     for i in range(1):
#             print('Optimisation trial {0}.'.format(i))
#             form = _form(form)
#             fopt, qopt = optimise_single(form, qmax=5, population=200, generations=200, printout=10, tol=0.01)
#             if fopt is not None and fopt < 2000:
#                 print('Optimisation found a result after {0} trials.'.format(i))
#                 break

#     form.to_json(file)
#     form_ = replicate(form,file_complete)

#     if plots:
#         plot_form(form, radius=0.05).show()
#         plot_form(form_, radius=0.05).show()

#     z_from_form(form_)
#     oveview_forces(form_)

#     if save:
#         form.to_json(file)
#         form_.to_json(file_complete)

#     return form, form_

# def find_independents(E):

#     _, m = E.shape
#     Etemp = E[:,[0]]
#     ind = []

#     for i in range(1,m):
#         Etest = hstack([Etemp,E[:,[i]]])
#         _ , ncol = Etest.shape
#         if matrix_rank(Etest) < ncol:
#             ind.append(i)
#         else:
#             Etemp = Etest

#     return ind

# def independents_exclude(E, outs):

#     _, m = E.shape
#     possible = list(set(range(m)) - set(outs))
#     Eouts = E[:,outs]
#     if matrix_rank(Eouts) == Eouts.shape[1]:
#         Etemp = Eouts
#         ind = []
#     else:
#         print('Warning, could not exclude all')
#         return find_independents(E)

#     for i in possible:
#         Etest = hstack([Etemp,E[:,[i]]])
#         _ , ncol = Etest.shape
#         if matrix_rank(Etest) < ncol:
#             ind.append(i)
#         else:
#             Etemp = Etest

#     return ind

# def independents_include(E, ins):

#     n, m = E.shape
#     if len(ins) > (m-n):
#         print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
#         ins = ins[:(m-n)]
#     not_in = list(set(range(m)) - set(ins))
#     Ein = E[:,ins]
#     while matrix_rank(Ein) < Ein.shape[1]:
#         print('Warning, edges are dependent among each other')
#         ins = ins[matrix_rank(Ein)]
#     if len(ins) == (m-n):
#         Enot_in = E[:,not_in]
#         if matrix_rank(Enot_in) == Enot_in.shape[1]:
#             return ins
#         else:
#             print('Warning, edges do not form an independent set')
#             return find_independents(E)
#     ind = ins
#     Etemp = E[:,[not_in[0]]]
#     Etemp.shape

#     for i in not_in[1:]:
#         Etest = hstack([Etemp,E[:,[i]]])
#         _ , ncol = Etest.shape
#         if matrix_rank(Etest) < ncol:
#             ind.append(i)
#         else:
#             Etemp = Etest

#     return ind

# def inds_incl_excl(E, ins, outs):

#     n, m = E.shape
#     if len(ins) > (m-n):
#         print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
#         ins = ins[:(m-n)]
#     not_in = list(set(range(m)) - set(ins))
#     possible = list(set(not_in)-set(outs))
#     Ein = E[:,ins]
#     Eouts = E[:,outs]
#     while matrix_rank(Ein) < Ein.shape[1]:
#         print('Warning, included edges are dependent among each other')
#         ins = ins[matrix_rank(Ein)]
#     if len(ins) == (m-n):
#         Enot_in = E[:,not_in]
#         if matrix_rank(Enot_in) == Enot_in.shape[1]:
#             return ins
#         else:
#             print('Warning, edges do not form an independent set')
#             return find_independents(E)
#     if matrix_rank(Eouts) == Eouts.shape[1]:
#         Etemp = Eouts
#         ind = ins
#     else:
#         print('Warning, could not exclude all')
#         ind = ins
#         Etemp = E[:,[not_in[0]]]
#         possible = not_in[1:]
    
#     for i in possible:
#         Etest = hstack([Etemp,E[:,[i]]])
#         _ , ncol = Etest.shape
#         if matrix_rank(Etest) < ncol:
#             ind.append(i)
#         else:
#             Etemp = Etest

#     return ind

# def initialize_problem(form, indset = None, printout = None):

#     # Mapping

#     k_i  = form.key_index()
#     i_k  = form.index_key()
#     i_uv = form.index_uv()
#     uv_i = form.uv_index()

#     # Vertices and edges

#     n     = form.number_of_vertices()
#     m     = form.number_of_edges()
#     fixed = [k_i[key] for key in form.fixed()]
#     rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
#     edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
#     sym   = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
#     free  = list(set(range(n)) - set(fixed) - set(rol))

#     # Constraints

#     lb_ind = []
#     ub_ind = []
#     lb = []
#     ub = []
#     for key, vertex in form.vertex.items():
#         if vertex.get('lb', None):
#             lb_ind.append(k_i[key])
#             lb.append(vertex['lb'])
#         if vertex.get('ub', None):
#             ub_ind.append(k_i[key])
#             ub.append(vertex['ub'])
#     lb = array(lb)
#     ub = array(ub)
#     lb.shape = (len(lb),1)
#     ub.shape = (len(ub),1)

#     # Co-ordinates and loads

#     xyz = zeros((n, 3))
#     x   = zeros((n, 1))
#     y   = zeros((n, 1))
#     z   = zeros((n, 1))
#     s   = zeros((n, 1))
#     px  = zeros((n, 1))
#     py  = zeros((n, 1))
#     pz  = zeros((n, 1))
#     s   = zeros((n, 1))
#     w   = zeros((n, 1))

#     for key, vertex in form.vertex.items():
#         i = k_i[key]
#         xyz[i, :] = form.vertex_coordinates(key)
#         x[i]  = vertex.get('x')
#         y[i]  = vertex.get('y')
#         px[i] = vertex.get('px', 0)
#         py[i] = vertex.get('py', 0)
#         pz[i] = vertex.get('pz', 0)
#         s[i]  = vertex.get('target', 0)
#         w[i] = vertex.get('weight', 1.0) # weight used in case of fiting...

#     Wfree = diags(w[free].flatten())
#     xy = xyz[:, :2]
#     z = xyz[:, 2]

#     # C and E matrices

#     C   = connectivity_matrix(edges, 'csr')
#     Ci  = C[:, free]
#     Cf  = C[:, fixed]
#     Ct = C.transpose()
#     Cit = Ci.transpose()
#     E   = equilibrium_matrix(C, xy, free, 'csr').toarray()
#     uvw = C.dot(xyz)
#     U   = uvw[:, 0]
#     V   = uvw[:, 1]

#     print('Equilibrium Matrix Shape: ', E.shape)
#     start_time = time.time()

#     # Independent and dependent branches

#     if indset:
#         ind = []
#         for u, v in form.edges():
#             if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
#                 ind.append(uv_i[(u, v)])
#     else:
#         _, s, _ = svd(E)
#         ind = independents_exclude(E, sym)
    
#     k   = len(ind)
#     dep = list(set(range(m)) - set(ind))
#     elapsed_time = time.time() - start_time
#     print('Found {0} independents'.format(k))
#     print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

#     for u, v in form.edges():
#         form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

#     if printout:
#         print('Form diagram has {0} (RREF) or {1} (SVD) independent branches '.format(len(ind), m - len(s)))
#         print('Shape Equilibrium Matrix: ', E.shape)
#         print('Found {0} independents'.format(k))

#     # Set-up

#     try:
#         t = form.attributes['offset']
#     except:
#         t = None

#     lh     = normrow(C.dot(xy))**2
#     Edinv  = -csr_matrix(pinv(E[:, dep]))
#     Ei     = E[:, ind]
#     p      = vstack([px[free], py[free]])
#     q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

#     tol = 0.01
#     opt_max = False
#     planar = False
#     tension = False
#     target = False
#     anchors = []
#     b = None


#     args   = (q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b)
    
#     return args










