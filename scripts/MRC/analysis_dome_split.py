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
from compas.colors import Color

from numpy import array

import os

xspan = yspan = [0, 10.0]
thk = 0.50
discretisation = [16, 20]
xc = yc = (xspan[1] - xspan[0])/2
radius = (xspan[1] + xspan[0])/2

dome = Shape.create_dome(thk=thk, radius=radius, t=1.0)

form = FormDiagram.create_circular_radial_form(discretisation=discretisation)

vector_supports = []
vectors_plot = []
base_plot = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x - xc > 0.1:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))
    if x - xc < -0.1:
        dXbi = [-1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

dXb = array(vector_supports)

constraints = ['funicular', 'envelope', 'reac_bounds']

analysis = Analysis.create_compl_energy_analysis(form,
                                                 dome,
                                                 printout=True,
                                                 support_displacement=dXb,
                                                 max_iter=2000,
                                                 solver='SLSQP',
                                                 starting_point='loadpath')

analysis.optimiser.set_constraints(constraints)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

analysis.run()

# folder = os.path.join('/Users/mricardo/compas_dev/me/compl_energy/dome/split', form.parameters['type']) + '/'
# os.makedirs(folder, exist_ok=True)
# title = dome.datashape['type'] + '_' + form.parameters['type'] + '_discr_' + str(discretisation)
# address = folder + title + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
# if analysis.optimiser.exitflag == 0:
#     print('Saving form to:', address)
#     form.to_json(address)

plotter = TNOPlotter(form)
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    plotter.draw_vector(vector=vector, base=base)
plotter.draw_form()
plotter.draw_supports()
plotter.draw_cracks()
plotter.show()

plotter = TNOPlotter(form)
plotter.settings['color.edges.independent'] = Color.blue()
plotter.draw_form_independents()
plotter.draw_supports()
plotter.show()

view = Viewer(form, shape=dome)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form()
view.draw_cracks()
view.draw_shape()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()
