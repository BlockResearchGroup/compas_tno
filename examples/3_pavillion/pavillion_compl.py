from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import apply_bounds_reactions
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape

from compas.geometry import Vector
from compas.geometry import Point

from compas.geometry import Scale

from numpy import array

import os

xspan = yspan = [0, 10.0]
thk = 0.50
discretisation = 10
pavillion_original = Shape.create_pavillionvault(thk=thk, t=5.0)
xc = yc = (xspan[1] - xspan[0])/2

pavillionexp = Shape.create_pavillionvault(thk=thk, t=5.0, expanded=True)

intra = pavillionexp.intrados
extra = pavillionexp.extrados
middle = pavillionexp.middle

scale = Scale.from_factors([1.0, 1.0, 3.5/5.0])  # change apex height to 3.5
intra.transform(scale)
extra.transform(scale)
middle.transform(scale)

pavillion = Shape.from_meshes(intra, extra, middle)

view = Viewer(shape=pavillion)
view.draw_shape()
view.show()

# form = FormDiagram.create_ortho_form(fix='all', discretisation=discretisation)
# form = FormDiagram.create_cross_form(fix='all', discretisation=discretisation)
form = FormDiagram.create_cross_with_diagonal(fix='all', discretisation=discretisation)

apply_bounds_reactions(form, pavillion_original)
apply_envelope_from_shape(form, pavillion)
apply_selfweight_from_shape(form, pavillion)

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if abs(x - xspan[1]) < 1e-3:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(dXbi[0], dXbi[2], 0.0))
        base_plot.append(Point(x, y, z - 0.2))

    vector_supports.append(dXbi)

dXb = array(vector_supports)

# constraints = ['funicular', 'envelope', 'reac_bounds']

# analysis = Analysis.create_compl_energy_analysis(form,
#                                                  pavillion,
#                                                  printout=True,
#                                                  support_displacement=dXb,
#                                                  max_iter=2000,
#                                                  solver='IPOPT',
#                                                  starting_point='loadpath')

# analysis.optimiser.set_constraints(constraints)
# # analysis.apply_selfweight()
# # analysis.apply_envelope()
# # analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()

# analysis.run()

# folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/pavillion/wall_open', form.parameters['type']) + '/'
# os.makedirs(folder, exist_ok=True)
# title = pavillion.datashape['type'] + '_' + form.parameters['type'] + '_discr_' + str(discretisation)
# address = folder + title + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
# if analysis.optimiser.exitflag == 0:
#     print('Saving form to:', address)
#     form.to_json(address)

# plotter = TNOPlotter(form)
# for i in range(len(vectors_plot)):
#     vector = vectors_plot[i]
#     base = base_plot[i]
#     plotter.draw_vector(vector=vector, base=base)
# plotter.draw_form()
# plotter.draw_supports()
# plotter.draw_cracks()
# plotter.show()

pavillion_view = Shape.create_pavillionvault(thk=thk, t=0.0, expanded=True)

intra = pavillion_view.intrados
extra = pavillion_view.extrados
middle = pavillion_view.middle

scale = Scale.from_factors([1.0, 1.0, 3.5/5.0])  # change apex height to 3.5
intra.transform(scale)
extra.transform(scale)
middle.transform(scale)

pavillion_nice = Shape.from_meshes(intra, extra, middle)

view = Viewer(form, shape=pavillion_nice)
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()

view = Viewer(form, shape=pavillion_nice)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = -180
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()
