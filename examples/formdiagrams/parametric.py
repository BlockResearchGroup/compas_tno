from matplotlib.pyplot import plot
from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.utilities import move_pattern_to_origin
from compas_tno.utilities import mesh_remove_two_valent_nodes
from compas.datastructures import network_join_edges
from compas_plotters import Plotter

for lambd in [1.0]:

    form : FormDiagram = FormDiagram.create_parametric_form(discretisation=10, lambd=lambd)
    path = '/Users/mricardo/compas_dev/me/pattern/parametric/' + 'form_lambd_' + str(lambd) + '.json'
    force = None

    form.uv_index()

    plotter = Plotter()
    art = plotter.add(form)
    art.draw_vertexlabels()
    plotter.show()

    print('edges', form.number_of_edges())
    print('vertices', form.number_of_vertices())
    print('faces', form.number_of_faces())

    for edge in form.edges():
        ec = form.edge_coordinates(*edge)
        print(edge, ec)

    lengths = [form.edge_length(u, v) for u, v in form.edges()]
    areas = [form.face_area(key) for key in form.faces()]

    print('max / min lengths:', max(lengths), min(lengths))
    print('max / min areas:', max(areas), min(areas))

    force = reciprocal_from_form(form, plot=True)

    # form.to_json(path)

    plotter = TNOPlotter(form, force=force)
    plotter.draw_form(scale_width=False)
    # plotter.draw_force()
    # plotter.formartist.draw_vertexlabels()
    plotter.draw_supports()
    plotter.show()
