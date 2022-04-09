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
thk = 0.5
discretisation = [8, 12]
discretisation_shape = [2*discretisation[0], 2*discretisation[1]]
lambd0 = 0.1

# Parameters Optimisations

obj = 'max_load'
solver = 'IPOPT'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb', 'lambdh']
features = ['fixed']
starting_point = 'loadpath'
make_video = True
autodiff = False

# Create shape/diagram

dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.0)

form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius, diagonal=True, partial_diagonal='right')
# form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)

# plotter = Plotter()
# plotter.fontsize = 6
# artist = plotter.add(form)
# artist.draw_vertexlabels()
# plotter.zoom_extents()
# plotter.show()

# form = FormDiagram.create_circular_spiral_form(discretisation=discretisation, radius=radius)

# Maximum load magnitude

# max_load_mult = 600.0
# n = form.number_of_vertices()
# pzv = zeros((n, 1))
# pzv[0] = -1.0
# # pzv[31] = -1.0

apply_selfweight_from_shape(form, dome)
apply_horizontal_multiplier(form, lambd=lambd0, direction='x')
n = form.number_of_vertices()
phv = zeros((2*n, 1))
for i, vertex in enumerate(form.vertices()):
    phv[i] = form.vertex_attribute(vertex, 'px')
    phv[i*2] = form.vertex_attribute(vertex, 'py')

# plotter = TNOPlotter(form)
# plotter.draw_form(scale_width=False)
# plotter.draw_supports()
# plotter.show()

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
optimiser.settings['save_force_diagram'] = False
optimiser.settings['autodiff'] = autodiff
optimiser.settings['max_iterations'] = 500

optimiser.settings['max_lambd'] = 5.0
optimiser.settings['load_direction'] = phv

# Create analysis

analysis = Analysis.from_elements(dome, form, optimiser)
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
view = Viewer(form, dome)
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

plotter = TNOPlotter(form, dome)
plotter.show_solution()

view = Viewer(form, dome)
view.draw_force()
view.show_solution()

z = [form.vertex_attribute(key, 'z') for key in form.vertices()]

print('end min, max z:', min(z), max(z))

if make_video:

    from compas_tno.viewers import animation_from_optimisation
    from compas_tno.algorithms import reciprocal_from_form
    import compas_tno

    DATA_XFORM = compas_tno.get('Xform.json')
    DATA_XFORCE = compas_tno.get('Xforce.json')

    force = reciprocal_from_form(form)

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
