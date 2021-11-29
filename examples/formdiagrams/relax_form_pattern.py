from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import reciprocal_from_form
from compas_plotters import MeshPlotter
from compas_tno.algorithms import apply_sag

data = {
    'type': 'cross_with_diagonal',  # 'cross_fd', 'ortho_fd'
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 10,
    'fix': 'corners',
}

for bf in [1.0, 2.0, 5.0, 10.0, 20.0]:

    form = FormDiagram.from_library(data)

    for key in form.vertices():
        form.vertex_attribute(key, 'pz', -1.0)

    apply_sag(form, boundary_force=bf)

    plotter = MeshPlotter(form, figsize=(8, 8))
    plotter.draw_edges()
    plotter.draw_vertices(keys=form.fixed(), facecolor='000000')
    plotter.show()

    force = reciprocal_from_form(form)
    # force = ForceDiagram.from_formdiagram(form)

    plotter = MeshPlotter(force, figsize=(8, 8))
    plotter.draw_edges()
    plotter.show()
