from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import apply_sag
from compas_tno.utilities.form import slide_pattern_inwards
from compas_tno.utilities.form import displacement_map_4parabolas
from compas_tno.utilities import form_parabolic_slide
from compas.colors import Color

path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/FormDiagram-crossbraced.json'
form = FormDiagram.from_json(path)

# plot = TNOPlotter(form)
# plot.draw_form(scale_width=False)
# plot.draw_supports()
# plot.show()

for delta in [0.0, 0.5, 0.714]:
    form2 = FormDiagram.create_cross_with_diagonal(discretisation=14)
    # apply_sag(form2)
    slide_pattern_inwards(form2, delta=delta)
    # displacement_map_4parabolas(form2)

    plot = TNOPlotter(form2)
    plot.draw_form(scale_width=False, color=Color.black())
    plot.draw_supports(color=Color.red(), size=8)
    plot.show()

# plot = TNOPlotter(form2)
# plot.draw_form(scale_width=False)
# plot.draw_mesh(form, color=Color.grey())
# plot.draw_supports()
# plot.show()
