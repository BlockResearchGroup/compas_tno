from ipopt import minimize_ipopt
# minimize_ipopt(fobj, x0, args = args, constraints = [fconstr])
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR DOME --------------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 16]

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'center': [5.0, 5.0],
    'radius': radius,
    't' : 10.0
}

dome = Shape.from_library(data_shape)
print('Dome created!')

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

form = form.initialise_tna(plot=False, kmax=50)
# plot_form(form).show()

# --------------------- 3.1. Analyse by hand here ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'

# optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
# optimiser.data['variables'] = ['ind', 'zb']
# optimiser.data['objective'] = 'min'

optimiser.data['constraints'] = ['funicular']
optimiser.data['variables'] = ['ind']
optimiser.data['objective'] = 'loadpath'


optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
print(optimiser.data)

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
# plot_independents(form).show()
analysis.run()


# from compas_tno.algorithms import set_up_nonlinear_optimisation
# analysis = set_up_nonlinear_optimisation(analysis)

# form = analysis.form
# optimiser = analysis.optimiser
# solver = optimiser.data['solver']
# fobj = optimiser.fobj
# fconstr = optimiser.fconstr
# args = optimiser.args
# q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints = args
# i_uv = form.index_uv()
# i_k = form.index_key()
# k_i = form.key_index()
# bounds = optimiser.bounds
# x0 = optimiser.x0
# g0 = optimiser.g0
# plot = optimiser.data['plot']

# lower = [lw[0] for lw in bounds]
# upper = [up[1] for up in bounds]
# cu = [10e20]*len(g0)
# cl = [0]*len(g0)

# problem_obj = wrapper_ipopt()
# problem_obj.fobj = fobj
# problem_obj.fconstr = fconstr
# problem_obj.args = args
# problem_obj.bounds = bounds
# problem_obj.x0 = x0

# nlp = ipopt.problem(
#         n=len(x0),
#         m=len(g0),
#         problem_obj=problem_obj,
#         lb=lower,
#         ub=upper,
#         cl=cl,
#         cu=cu
#         )

# xopt, info = nlp.solve(x0)

# print(xopt)
# print(info)


# Find independent edges
analysis.run()
# plot_form(form, show_q=False).show()
# view_thrust(form).show()

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form.to_json(file_adress)
# form.from_json(file_adress)

# # --------------------- 4.1 Create Minimisation Optimiser ---------------------

# optimiser = Optimiser()
# optimiser.data['library'] = 'MMA'
# optimiser.data['solver'] = 'MMA'
# optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
# optimiser.data['variables'] = ['ind', 'zb']
# optimiser.data['objective'] = 'min'
# optimiser.data['printout'] = True
# optimiser.data['solver_options']['derivatives'] = 'DF_reduced'  # 'DF_brute' 'DF_reduced' and 'analytical' in process.
# optimiser.data['plot'] = True
# optimiser.data['find_inds'] = True
# optimiser.data['qmax'] = 50.0
# print(optimiser.data)

# # --------------------- 4.2 Create Minimisation Optimiser ---------------------

# analysis = Analysis.from_elements(dome, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser() # Find independent edges
# analysis.run()

# form = analysis.form
plot_form(form, show_q=False).show()

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form.to_json(file_adress)

view_thrust(form).show()
