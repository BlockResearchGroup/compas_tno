from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter

path = '/Users/mricardo/compas_dev/me/compl_energy/dome/split/radial_fd/dome_radial_fd_discr_[16, 20]_Ecomp-linear_thk_50.0.json'

form = FormDiagram.from_json(path)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.draw_supports()
plotter.draw_reactions()
plotter.draw_vertexlabels()
plotter.show()
