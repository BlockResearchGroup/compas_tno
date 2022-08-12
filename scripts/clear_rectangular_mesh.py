from compas.datastructures import Mesh
from compas.geometry import Translation
from compas.geometry import Transformation
from compas.geometry import distance_point_point_xy
from compas.geometry import Vector
from compas_plotters import Plotter
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter

from compas_view2.app import App

# jsonpath = '/Users/mricardo/compas_dev/me/freeform/meshes_rectangular/'
# intra_file = jsonpath + '1-intrados-mesh.json'
# extra_file = jsonpath + '1-extrados-mesh.json'

# intrados = NurbsSurface.from_json(intra_file)
# extrados = NurbsSurface.from_json(extra_file)

# intrados = Mesh.from_json(intra_file)
# extrados = Mesh.from_json(extra_file)

# print(intrados)

# app = App(show_grid=False)
# app.add(intrados)
# app.add(extrados)
# app.show()


# load and interpret meshes

for prob in ['E']:  # ['A', 'B', 'C', 'D', 'E']

    mesh_file = '/Users/mricardo/compas_dev/me/freeform/meshes_square/mesh-' + prob + '.json'
    form_file = '/Users/mricardo/compas_dev/me/freeform/meshes_square/form-' + prob + '.json'

    mesh_file = '/Users/mricardo/compas_dev/compas_tno/data/CISM/test-CROSSFAN.json'
    form_file = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CROSSFAN.json'

    mesh = Mesh.from_json(mesh_file)
    print('nedges', mesh.number_of_edges())

    bbox = mesh.bounding_box_xy()

    xmin = min([bbox[i][0] for i in range(len(bbox))])
    ymin = min([bbox[i][1] for i in range(len(bbox))])
    xmax = max([bbox[i][0] for i in range(len(bbox))])
    ymax = max([bbox[i][1] for i in range(len(bbox))])

    dx = (xmax - xmin)
    dy = (ymax - ymin)

    vec = Vector(-xmin, -ymin, 0)
    T = Translation.from_vector(vec)
    mesh = mesh.transformed(T)

    bbox = mesh.bounding_box_xy()
    # print('New bbox:', bbox)

    plotter = Plotter()
    artist = plotter.add(mesh)
    artist.draw_vertexlabels()
    plotter.show()

    form = FormDiagram.from_mesh(mesh)

    for vertex in form.vertices():
        coord = form.vertex_coordinates(vertex)
        dist_corners = min([
            distance_point_point_xy(coord, [0.0, 0.0]),
            distance_point_point_xy(coord, [0.0, dy]),
            distance_point_point_xy(coord, [dx, 0.0]),
            distance_point_point_xy(coord, [dx, dy])
        ])
        if dist_corners < 1e-3:
            form.vertex_attribute(vertex, 'is_fixed', True)

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_supports()
    plotter.show()

    plotter = Plotter()
    artist = plotter.add(mesh)
    artist.draw_vertexlabels()
    plotter.show()

    dels = []
    for face in form.faces():
        area = form.face_area(face)
        if area > 10:
            dels.append(face)
            print(face, area)

    for delete in dels:
        form.delete_face(delete)

    plotter = Plotter()
    artist = plotter.add(mesh)
    artist.draw_vertexlabels()
    plotter.show()

    form.to_json(form_file)
    print('Form Saved to:', form_file)



