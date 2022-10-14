from compas_tno import analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer

discr = [20, 16]
# discr = [12, 12]

# path = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_compression_neg.json'
# form = FormDiagram.from_json(path)

form = FormDiagram.create_circular_radial_form(discretisation=discr)
shape = Shape.create_dome(thk=0.50, t=0.0)
shape.ro = 20.0

analysis = Analysis.create_minthk_analysis(form,
                                           shape,
                                           printout=True,
                                           solver='IPOPT',
                                           starting_point='loadpath',
                                           plot=True
                                           )
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.optimiser.set_features(['fixed'])

analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_bounds_on_q(qmax=1e-6)
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

# pz = {key: form.vertex_attribute(key, 'b') for key in form.vertices()}
# plotter = TNOPlotter(form, shape)
# plotter.draw_form(scale_width=False)
# plotter.draw_vertexlabels(text=pz)
# plotter.show()

analysis.run()

M = analysis.optimiser.M

# pz = {key: form.vertex_attribute(key, 'b') for key in form.vertices()}
# plotter = TNOPlotter(form, shape)
# plotter.draw_form(scale_width=False)
# plotter.draw_vertexlabels(text=pz)
# plotter.show()

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
