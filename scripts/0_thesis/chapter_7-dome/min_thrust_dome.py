from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.solvers.solver_MMA import analytical_f_g
from compas_tno.viewers import Viewer

rho = 1.0
thk = 0.5
discr = [20, 16]
# discr = [12, 12]

# path = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_compression_neg.json'
# form = FormDiagram.from_json(path)

form = FormDiagram.create_circular_radial_form(discretisation=discr)
shape = Shape.create_dome(thk=thk, t=1.0)
shape.ro = rho

# analysis = Analysis.create_minthrust_analysis(form, shape, printout=True, solver='IPOPT')
analysis = Analysis.create_maxthrust_analysis(form, shape, printout=True, solver='IPOPT')
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.optimiser.set_features(['fixed'])

analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_bounds_on_q(qmax=1e-6)
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

fopt = analysis.optimiser.fopt
W = form.lumped_swt()

obj = analysis.optimiser.settings['objective']

print('fopt/W:', fopt/W)

# pz = {key: round(form.vertex_attribute(key, 'pz'), 1) for key in form.vertices()}
# plotter = TNOPlotter(form, shape)
# plotter.draw_form(scale_width=False)
# plotter.draw_vertexlabels(text=pz)
# plotter.show()

# plotter = TNOPlotter(form, shape)
# plotter.draw_form()
# plotter.draw_cracks()
# plotter.draw_supports()
# plotter.show()

shape = Shape.create_dome()
view = Viewer(form, shape=shape)
view.scale_edge_thickness(10/rho)
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

folder = '/Users/mricardo/compas_dev/me/minmax/dome/thesis/'
# title = 'dome_{}_thk_{}_{}.json'.format(discr, thk, analysis.optimiser.settings['objective'])
# save = folder + title
# form.to_json(save)
# print('Saved to:', save)

title_raw = 'dome_{}_thk_{}_{}_raw.json'.format(discr, thk, analysis.optimiser.settings['objective'])
save_raw = folder + title_raw
view.to_json(save_raw)
