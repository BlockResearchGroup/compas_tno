from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution
from compas_tno.problems import initialise_problem
from compas_tno.algorithms import zlq_from_q
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import f_min_thrust_general
from compas_tno.problems import constr_wrapper_general
from compas_tno.problems import sensitivities_wrapper_general
from compas_tno.problems import gradient_fmin_general

from compas_plotters import MeshPlotter

from compas_tno.algorithms import apply_sag

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis

from numpy import zeros
from numpy import ones

from numpy import hstack
from numpy import vstack


span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.5
discretisation_shape = 4 * discretisation

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)
# plot_form(form, show_q=False).show()

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    # 'hc': hc,
    # 'hm': None,
    # 'he': None,
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 0.0,
}

vault = Shape.from_library(data_shape)

# -----------
# ----------- SYMMETRY -------------


from compas_tno.problems.setup import set_symmetry_constraint

# Horizontal Symmetry
# a = [0.0, 5.0, 0.0]
# b = [10.0, 5.0, 0.0]

# Vertical Symmetry
# a = [5.0, 0.0, 0.0]
# b = [5.0, 10.0, 0.0]

# Diagonal Symmetry
# a = [5.0, 5.0, 0.0]
# b = [0.0, 0.0, 0.0]

# form.apply_symmetry(axis_symmetry=[a, b])

# ----------------------------------------------------------
# ------------------------ SYMMETRY ------------------------
# ----------------------------------------------------------

form.apply_symmetry()
Esym = form.build_symmetry_transformation(printout=True)
Asym = form.build_symmetry_matrix(printout=False)
Asup = form.build_symmetry_matrix_supports(printout=False)

# Aglue = form.assemble_symmetry_matrix()

print(Asym.shape)
print(Esym.shape)

from compas_tno.algorithms import find_independents

ind = find_independents(Asym)
print(ind)
print(len(ind))

from compas_tno.plotters import plot_independents
from compas_tno.plotters import plot_sym_inds

# plot_sym_inds(form).show()

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

# Apply Selfweight and Envelope

form.envelope_from_shape(vault)
form.selfweight_from_shape(vault)

form_base = form.copy()

# apply_sag(form)

# form.initialise_loadpath()

form.envelope_on_x_y(c=0.50)

for key in form.vertices_where({'is_fixed': True}):
    form.vertex_attribute(key, 'z', 1.00)

# plot_form(form, show_q=False).show()
# view_solution(form).show()

# ------------------------------------------------------------
# ----------------------- DOING BY HAND ----------------------
# ------------------------------------------------------------

# Extract the matrices

# form.edges_attribute('q', 25.0)
# for edge in form.edges_on_boundary():
#     form.edge_attribute(edge, 'q', 50.0)
# # args = initialise_problem(form, find_inds=False)

from compas_tno.problems import Problem
M = Problem.from_formdiagram(form)

# M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cf)

# # Variables are q [n x 1] and Xb
# variables = vstack([M.q, M.X[M.fixed].flatten('F').reshape((-1, 1))])

# hor = f_min_thrust_general(variables, M)
# print('Horizontal reaction:', hor)

# constr = constr_wrapper_general(variables, M)
# # print(constr)
# print('Constraint limits:', min(constr), max(constr))
# print('Shape constraint vector:', constr.shape)
# print('Length of Variables:', len(variables))

# deriv = sensitivities_wrapper_general(variables, M)
# # print('deriv limits:', min(deriv), max(deriv))
# print('deriv shape:', deriv.shape)

# grad = gradient_fmin_general(variables, M)
# print('grad shape:', grad.shape)

# i = 0
# for key in form.vertices():
#     form.vertex_attribute(key, 'x', M.X[i, 0])
#     form.vertex_attribute(key, 'y', M.X[i, 1])
#     form.vertex_attribute(key, 'z', M.X[i, 2])
#     # zfree.append(M.X[i, 2])
#     i = i + 1

E = M.E
print(E.shape)

from numpy import vstack
import time

i_uv = form.index_uv()
uv_i = form.uv_index()

# --------------------------
print('\n')

start_time = time.time()

Ass = vstack([E, Asym])
print(Ass.shape)

ind = find_independents(Ass)
print(ind)
print(len(ind))

elapsed_time = time.time() - start_time
print('Elapsed Time blocks idea:', elapsed_time)
print('Number of independentes:', len(ind))

for i in ind:
    u, v = i_uv[i]
    print(i)
    form.edge_attribute((u, v), 'is_ind', True)
plot_independents(form).show()
for u, v in form.edges():
    form.edge_attribute((u, v), 'is_ind', False)

# --------------------------
print('\n')

start_time = time.time()

ind = find_independents(Asym)
print(ind)
print(len(ind))

elapsed_time = time.time() - start_time
print('Elapsed Time on symmetric matrix:', elapsed_time)
print('Number of independentes:', len(ind))

for i in ind:
    u, v = i_uv[i]
    form.edge_attribute((u, v), 'is_ind', True)
plot_independents(form).show()
for u, v in form.edges():
    form.edge_attribute((u, v), 'is_ind', False)

# --------------------------
print('\n')

start_time = time.time()
ind = find_independents(E)
print(ind)
print(len(ind))

elapsed_time = time.time() - start_time
print('Elapsed Time normal:', elapsed_time)
print('Number of independentes:', len(ind))

for i in ind:
    u, v = i_uv[i]
    form.edge_attribute((u, v), 'is_ind', True)
plot_independents(form).show()
for u, v in form.edges():
    form.edge_attribute((u, v), 'is_ind', False)

# --------------------------
print('\n')

start_time = time.time()

Ass = E.dot(Esym)
print(Ass.shape)

ind_sym = find_independents(Ass)
print(ind_sym)
print(len(ind_sym))

ind = []
for i in ind_sym:
    for u, v in form.edges():
        if form.edge_attribute((u, v), 'sym_key') == i:
            ind.append(uv_i[(u, v)])
            break

print('Make the equivalent of the indices in the network:')
print(ind)
print(len(ind))

elapsed_time = time.time() - start_time
print('Elapsed Time new Sym Matrix:', elapsed_time)
print('Number of independentes:', len(ind))

for i in ind:
    u, v = i_uv[i]
    form.edge_attribute((u, v), 'is_ind', True)
plot_independents(form).show()
for u, v in form.edges():
    form.edge_attribute((u, v), 'is_ind', False)

# ------------------

length_dep = Esym.shape[1]
dep = list(set(range(length_dep)) - set(ind))
print(dep)

from scipy.sparse import csr_matrix
from numpy.linalg import pinv

Edinv = -csr_matrix(pinv(Ass[:, dep]))
Ei = E[:, ind]
print(Edinv.shape)
print(Ei.shape)

Ehor = Edinv.dot(Ei)
print(Ehor.shape)


# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
optimiser.data['library'] = 'SLSQP'
optimiser.data['solver'] = 'SLSQP'
# optimiser.data['library'] = 'IPOPT'
# optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['q', 'Xb']
# optimiser.data['variables'] = ['q']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = False
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = True
optimiser.data['jacobian'] = True
optimiser.data['derivative_test'] = True

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# print(max(zfree), min(zfree))

for key in form.vertices_where({'is_fixed': True}):
    print(form.vertex_coordinates(key))

plot_superimposed_diagrams(form, form_base).show()

# Viewing
plot_form(form, show_q=False).show()
view_solution(form).show()
