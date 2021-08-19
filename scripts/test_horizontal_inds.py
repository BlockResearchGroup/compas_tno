
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
from compas_tno.shapes import Shape
from compas_tno.viewers import view_thrust
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter



# --------------------- 1. Parameters ---------------------

L = 10.0
discr = [4,4]
target_x  = L/2
target_y = L
load_mult_x = 1.0
load_mult_y = 1.0

# --------------------- 1. Create Form ---------------------

data = {
    'type': 'cross_fd',
    'xy_span': [[0.0, L], [0.0, L]],
    'discretisation': discr,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
# form = form.delete_boundary_edges()
# plot_form(form, show_q=False).show()

# ---------- 1. Create Pavillion Shape - to Assign loads ----------

data_shape = {
    'type': 'pavillionvault',
    'thk': 0.20,
    'discretisation': discr,
    'xy_span': [[0.0, L], [0.0, L]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()
print('Total SWT:', swt)

# --------------------- 3. Create Convex Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'MATLAB'
optimiser.settings['solver'] = 'SDPT3'
optimiser.settings['constraints'] = ['funicular']
optimiser.settings['variables'] = ['ind']
optimiser.settings['objective'] = 'loadpath'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 150.0

# -------------- 4. Create Analysis Model and Run Convex Opt --------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()

for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    if abs(x - target_x) < 10e-4 and abs(y - target_y) < 10e-4:
        target_key = key
try:
    print('Key(s) Loaded: ', target_key)
except BaseException:
    print('Did not get any key!')

target_x  = L/2
target_y = L/2

for key in form.vertices():
    x, y, z = form.vertex_coordinates(key)
    if abs(x - target_x) < 10e-4 and abs(y - target_y) < 10e-4:
        target_key2 = key
try:
    print('Key(s) Loaded: ', target_key2)
except BaseException:
    print('Did not get any key!')

analysis.apply_pointed_load(target_key, load_mult_x, component='px')
analysis.apply_pointed_load(target_key, load_mult_y, component='py')
analysis.apply_pointed_load(target_key2, load_mult_x, component='px')
analysis.apply_pointed_load(target_key2, load_mult_y, component='py')

plotter = MeshPlotter(form, figsize=(10,10))
plotter.draw_edges(keys=[key for key in form.edges_where({'_is_edge': True})])
plotter.draw_vertices(radius=0.10)
plotter.draw_vertices(keys=[target_key], text='px= {0:.1} py={1:.1}'.format(form.vertex_attribute(target_key, 'px'), form.vertex_attribute(target_key, 'py')), radius=0.10, facecolor='FF0000')
plotter.draw_vertices(keys=[target_key2], text='px= {0:.1} py={1:.1}'.format(form.vertex_attribute(target_key2, 'px'), form.vertex_attribute(target_key2, 'py')), radius=0.10, facecolor='FF0000')
plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.10, facecolor='000000')
plotter.show()  # Simple plot of the diagram


# -------------- 5. Set Up Optimiser -------------------

analysis.set_up_optimiser()
plot_independents(form).show()
analysis.run()

q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i = optimiser.args
# for u, v in form.edges_where('_is_edge': True):

import numpy
from numpy.linalg import matrix_rank
from numpy import vstack
from numpy import hstack
from numpy import zeros

numpy.set_printoptions(threshold=10000, linewidth=1000)
k = E.shape[0]
for i in range(k):
    new_vec = zeros((k,1))
    new_vec[i] = 1.0
    rank = matrix_rank(hstack([E,new_vec]))
    if rank
    print(i)
    print(matrix_rank(hstack([E,new_vec])))

plot_form(form).show()
view_thrust(form).show()
