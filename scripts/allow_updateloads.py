from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer


# Geometry parameters

radius = 5.0
thk = 0.5
discretisation = [8, 10]
discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

# Parameters Optimisations

obj = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']  # , 'update-loads']
starting_point = 'tna'
make_video = True

# Create shape/diagram

dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=1.0)

form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)

# Create optimiser

optimiser = Optimiser()
optimiser.settings['objective'] = obj
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['starting_point'] = starting_point
optimiser.settings['derivative_test'] = True
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['save_iterations'] = make_video

# Create analysis

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

analysis.run()

view = Viewer(form, dome)
view.draw_thrust()
view.draw_shape()
view.draw_force()
view.draw_cracks()
view.draw_reactions()
view.show()


# # # ---- MAKE THE VIDEO ----

# if make_video:

#     from compas_tno.viewers import animation_from_optimisation
#     from compas_tno.algorithms import reciprocal_from_form
#     import compas_tno

#     DATA_XFORM = compas_tno.get('Xform.json')
#     DATA_XFORCE = compas_tno.get('Xforce.json')

#     force = reciprocal_from_form(form)

#     animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, interval=150)
