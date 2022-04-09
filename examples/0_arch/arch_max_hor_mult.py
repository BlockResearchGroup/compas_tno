from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter

from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_selfweight_from_shape

from compas_plotters import Plotter
from numpy import zeros
import os

# Geometry parameters

radius = 5.0
thk = 0.1
span = 1.0
discretisation = 20
lambd0 = 0.1

H = span/2
L = span

# Parameters Optimisations

obj = 'max_load'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb', 'lambdh']
features = ['fixed']
starting_point = 'loadpath'
make_video = True
show_force_diagram = False
autodiff = False

arch = Shape.create_arch(H=H, L=L, thk=thk, t=0.5, discretisation=discretisation)

form = FormDiagram.create_arch(H=span/2, L=span, discretisation=discretisation)

apply_selfweight_from_shape(form, arch)
apply_horizontal_multiplier(form, lambd=lambd0, direction='x')
n = form.number_of_vertices()
phv = zeros((2*n, 1))
for i, vertex in enumerate(form.vertices()):
    phv[i] = form.vertex_attribute(vertex, 'px')
    phv[i*2] = form.vertex_attribute(vertex, 'py')

# Create optimiser

optimiser = Optimiser()
optimiser.settings['objective'] = obj
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['starting_point'] = starting_point
optimiser.settings['derivative_test'] = False
optimiser.settings['printout'] = True
optimiser.settings['plot'] = True
optimiser.settings['save_iterations'] = make_video
optimiser.settings['save_force_diagram'] = show_force_diagram
optimiser.settings['autodiff'] = autodiff
optimiser.settings['max_iterations'] = 500

optimiser.settings['max_lambd'] = 10.0
optimiser.settings['load_direction'] = phv

# Create analysis

analysis = Analysis.from_elements(arch, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()

pz0 = form.vertex_attribute(0, 'pz')

pzt = 0
for key in form.vertices():
    pz = form.vertex_attribute(key, 'pz')
    pzt += pz

print('Total load of:', pzt)

analysis.set_up_optimiser()

# to view starting point
view = Viewer(form, arch)
view.draw_thrust()
view.draw_shape()
view.show()

analysis.run()

fopt = analysis.optimiser.fopt
exitflag = analysis.optimiser.exitflag

lambd = fopt * lambd0

print('Exitflag is:', exitflag)

pc = lambd/pzt

print('Optimum horizontal multiplier is:', round(pc*100, 3), '%')

# folder = os.path.join('/Users/mricardo/compas_dev/me/max_load/dome/apex/', 'dome', 'radial_fd')
# os.makedirs(folder, exist_ok=True)
# title = 'dome' + '_' + 'radial_fd' + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective'] + '_thk_' + str(100*thk) + '_pct_stw_' + str(pc) + '.json'
# save_form = os.path.join(folder, title)
# form.to_json(save_form)

# print('Solution Saved at:', save_form)

plotter = TNOPlotter(form, arch)
plotter.settings['scale.reactions'] = 0.15
plotter.draw_form_xz()
plotter.draw_shape_xz()
plotter.draw_cracks()
plotter.draw_reactions()
plotter.show()

view = Viewer(form, arch)
# view.draw_force()
view.show_solution()

if make_video:

    from compas_tno.viewers import animation_from_optimisation
    from compas_tno.algorithms import reciprocal_from_form
    import compas_tno

    DATA_XFORM = compas_tno.get('Xform.json')
    DATA_XFORCE = compas_tno.get('Xforce.json')

    if show_force_diagram:
        force = reciprocal_from_form(form)
    else:
        DATA_XFORCE = None
        force = None

    jump = 1
    import json
    with open(DATA_XFORM, mode='r', encoding='utf-8') as f:
        Xform = json.load(f)
    if len(Xform) > 100:
        jump = int(len(Xform)/50)
    if len(Xform) > 300:
        jump = int(len(Xform)/100)
    if len(Xform) > 500:
        jump = int(len(Xform)/200)

    animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, interval=150, jump_each=jump)
