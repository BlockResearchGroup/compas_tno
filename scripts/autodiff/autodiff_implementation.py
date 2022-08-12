from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

# Parameters Geometries

type_formdiagram = 'radial_fd'
type_structure = 'dome'
thk = 0.5
discretisation = [8, 10]
discretisation_shape = [2 * discretisation[0], 2 * discretisation[1]]

# Parameters Optimisations

obj = 'min'
solver = 'SLSQP'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'loadpath'

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': 5.0,
    'discretisation': discretisation
}

form = FormDiagram.from_library(data_diagram)

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'center': [5.0, 5.0],
    'radius': 5.0,
    't': 0.0,
}

dome = Shape.from_library(data_shape)

# view = Viewer(form, dome)
# view.draw_thrust()
# view.draw_shape()
# view.show()

opt = Optimiser()
opt.settings['objective'] = obj
opt.settings['solver'] = solver
opt.settings['constraints'] = constraints
opt.settings['variables'] = variables
opt.settings['features'] = features
opt.settings['starting_point'] = starting_point
opt.settings['printout'] = True

print(dome.datashape)

analysis = Analysis.from_elements(dome, form, opt)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

# ------------

M = opt.M
x0 = opt.x0


from compas_tno.problems import f_min_thrust
from compas_tno.problems import gradient_fmin

print('-'* 10, 'FROM ANALYTIC', '-'* 10)
pt1 = f_min_thrust(x0, M)
a1 = gradient_fmin(x0, M)
print(a1)
print(pt1)

# make non sparse

M.C = M.C.toarray()
M.Ci = M.Ci.toarray()
M.Cit = M.Cit.toarray()
M.Cb = M.Cb.toarray()

from compas_tno.autodiff.jax_grad import xyz_from_q_jax
from compas_tno.autodiff.jax_grad import f_min_thrust_jax
from jax import grad

# grad_xyz = grad(xyz_from_q_jax)
grad_min = grad(f_min_thrust_jax)

print('-'* 10, 'FROM JAX', '-'* 10)
pt2 = f_min_thrust_jax(x0, M)
a2 = grad_min(x0, M)
print(a2)
print(pt2)
# b = grad_xyz(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)


print(a1.shape)
print(a2.shape)

diff = abs(a1.reshape(-1, 1) - a2)
print(max(diff))

# print(a)

# analysis.run()

# view = Viewer(form, dome)
# view.draw_thrust()
# view.draw_shape()
# view.draw_cracks()
# view.show()
