from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer

path = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_review.json'
form = FormDiagram.from_json(path)

plt = TNOPlotter(form)
plt.draw_form()
plt.draw_supports()
plt.draw_reactions()
plt.show()

view = Viewer(form)
view.draw_form()
view.draw_shape()
view.draw_reactions(emerging_reactions=True)
view.show()
