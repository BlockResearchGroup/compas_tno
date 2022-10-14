## Problem of minimum thickness of the dome with discr = [20, 16]
## Analysed with IPOPT this results in a solution dependent on rho
## Analysed with SLSQP this can not be solved

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter

discr = [20, 16]
rho = 20.0
thk = 0.50

form = FormDiagram.create_circular_radial_form(discretisation=discr)
shape = Shape.create_dome(thk=thk, t=1.0)
shape.ro = rho

# for key in form.vertices_where({'is_fixed': True}):
#     form.vertex_attribute(key, 'z', 0.1)

analysis = Analysis.create_minthk_analysis(form,
                                           shape,
                                           printout=True,
                                           solver='IPOPT',
                                           starting_point='loadpath',
                                           plot=False,
                                           max_iter=1500,
                                           )
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.optimiser.set_features(['fixed'])
analysis.apply_selfweight()
analysis.apply_envelope()
# analysis.optimiser.set_additional_options(solver_convex='MATLAB')
analysis.optimiser.set_additional_options(jacobian=True)
# analysis.apply_bounds_on_q(qmax=1e-6)
analysis.apply_reaction_bounds()
# analysis.set_optimiser_options(nlp_scaling_method='gradient-based')
analysis.set_up_optimiser()

pz = {key: round(form.vertex_attribute(key, 'lb'), 6) for key in form.vertices_where({'is_fixed': True})}
plotter = TNOPlotter(form, shape)
plotter.draw_form(scale_width=False)
plotter.draw_vertexlabels(text=pz)
plotter.show()

analysis.run()

M = analysis.optimiser.M

pz = {key: round(form.vertex_attribute(key, 'lb'), 6) for key in form.vertices_where({'is_fixed': True})}
plotter = TNOPlotter(form, shape)
plotter.draw_form(scale_width=False)
plotter.draw_vertexlabels(text=pz)
plotter.show()

# plotter = TNOPlotter(form, shape)
# plotter.draw_form()
# plotter.draw_cracks()
# plotter.draw_supports()
# plotter.show()

view = Viewer(form, shape=shape)
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
