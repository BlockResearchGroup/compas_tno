## Problem of minimum thickness of the dome with discr = [20, 16]
## Analysed with IPOPT this results in a solution dependent on rho
## Analysed with SLSQP this can not be solved

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
import math

discr = 14
rho = 1.0
thk = 0.50

solver = 'IPOPT'

span = 10.0
spr_angle = 30.0
alpha = 1/math.cos(math.radians(spr_angle))
L = span * alpha
Ldiff = L - span
xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]

form = FormDiagram.create_cross_form(discretisation=discr)
shape = Shape.create_crossvault(xy_span=xyspan_shape, thk=thk)
shape.ro = rho

# for key in form.vertices_where({'is_fixed': True}):
#     form.vertex_attribute(key, 'z', 0.1)

########### ------------------------
########### min thk problem
########### ------------------------

analysis = Analysis.create_minthk_analysis(form,
                                           shape,
                                           printout=True,
                                           solver=solver,
                                           starting_point='loadpath',
                                           plot=False
                                           )
analysis.optimiser.set_constraints(['funicular', 'envelope'])
analysis.optimiser.set_features(['fixed'])
analysis.apply_selfweight()
analysis.apply_envelope()
# analysis.apply_bounds_on_q(qmax=1e-6)
# analysis.apply_reaction_bounds()
# analysis.set_optimiser_options(nlp_scaling_method='gradient-based')
analysis.set_up_optimiser()

analysis.run()

shape_plot = Shape.from_formdiagram_and_attributes(form)
view = Viewer(form, shape=shape_plot)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()
view.show()

########### ------------------------
########### min thrust problem
########### ------------------------


analysis = Analysis.create_minthrust_analysis(form,
                                              shape,
                                              printout=True,
                                              solver=solver,
                                              starting_point='loadpath',
                                              plot=False
                                              )
analysis.optimiser.set_constraints(['funicular', 'envelope'])
analysis.optimiser.set_features(['fixed'])
analysis.apply_selfweight()
analysis.apply_envelope()
# analysis.apply_bounds_on_q(qmax=1e-6)
# analysis.apply_reaction_bounds()
# analysis.set_optimiser_options(nlp_scaling_method='gradient-based')
analysis.set_up_optimiser()

analysis.run()

print('thrust over weight', form.thrust()/form.lumped_swt())

shape_plot = Shape.from_formdiagram_and_attributes(form)
view = Viewer(form, shape=shape_plot)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()
view.show()

########### ------------------------
########### max hor load problem
########### ------------------------

from numpy import zeros
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import slide_diagram

apply_horizontal_multiplier(form, lambd=0.1, direction='x')
slide_diagram(form, delta=-0.5, tappered=True)

n = form.number_of_vertices()
load_direction = zeros((2*n, 1))
for i, vertex in enumerate(form.vertices()):
    load_direction[i] = form.vertex_attribute(vertex, 'px')
    load_direction[i*2] = form.vertex_attribute(vertex, 'py')

problem = Analysis.create_max_load_analysis(form, shape,
                                            horizontal=True,
                                            solver=solver,
                                            load_direction=load_direction,
                                            max_lambd=5.0,
                                            plot=True,
                                            printout=True,
                                            max_iter=5000
                                            )
problem.apply_envelope()
problem.optimiser.set_features(['fixed', 'sym'])
axis_sym = [[[0, 5.0], [10, 5.0]]]
problem.optimiser.set_axis_symmetry(axis_sym)
problem.set_up_optimiser()
problem.run()

shape_plot = Shape.from_formdiagram_and_attributes(form)
view = Viewer(form, shape=shape_plot)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 60
view.draw_form(cull_negative=True)
view.draw_cracks(cull_negative=True)
view.draw_reactions(extend_reactions=True)
view.draw_shape()
view.show()
