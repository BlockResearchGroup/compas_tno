from compas.datastructures import Mesh
from compas.datastructures import mesh_dual

from compas.colors import Color

from compas_view2.app import App

from compas_assembly.datastructures import Assembly

from compas.geometry import add_vectors, scale_vector


def extended_dual(form, cls=None):
    """Create the extended dual of the mesh, which is the centroid dual added with the faces in the boundary.

    Parameters
    ----------
    form : FormDiagram
        Form Diagram as the primal
    cls : Class, optional
        The class which should return the dual, by default None

    Returns
    -------
    dual
        Extended dual mesh in the class required
    """

    dual = mesh_dual(form, cls)
    dual.flip_cycles()

    edge_vertex = {}

    boundary_edges = [edge for boundary in form.edges_on_boundaries() for edge in boundary]

    for u, v in boundary_edges:
        x, y, z = form.edge_midpoint(u, v)
        edge_vertex[u, v] = edge_vertex[v, u] = dual.add_vertex(x=x, y=y, z=z)

    boundary_vertices = [vertex for boundary in form.vertices_on_boundaries() for vertex in boundary]
    dual_vertices = list(dual.vertices())

    for vertex in boundary_vertices:
        x, y, z = form.vertex_coordinates(vertex)
        vertices = [dual.add_vertex(x=x, y=y, z=z)]
        nbrs = form.vertex_neighbors(vertex, ordered=True)[::-1]

        for nbr in nbrs:
            if form.is_edge_on_boundary(vertex, nbr):
                vertices.append(edge_vertex[vertex, nbr])

            face = form.halfedge[vertex][nbr]
            if face is not None:
                if face in dual_vertices:
                    vertices.append(face)

        dual.add_face(vertices)

    return dual


def offset_dual(dual, thk):
    idos = dual.copy()
    edos = dual.copy()

    assembly = Assembly()

    for vertex in dual.vertices():
        point = dual.vertex_coordinates(vertex)
        normal = dual.vertex_normal(vertex)

        idos.vertex_attributes(vertex, 'xyz', add_vectors(point, scale_vector(normal, -0.5 * thk)))
        edos.vertex_attributes(vertex, 'xyz', add_vectors(point, scale_vector(normal, 0.5 * thk)))

    for face in idos.faces():
        bottom = idos.face_coordinates(face)
        top = edos.face_coordinates(face)

        f = len(bottom)

        faces = [
            list(range(f)),
            list(range(f + f - 1, f - 1, -1))]

        for i in range(f - 1):
            faces.append([i, i + f, i + f + 1, i + 1])
        faces.append([f - 1, f + f - 1, f, 0])

        block = Mesh.from_vertices_and_faces(bottom + top, faces)
        assembly.add_block(block)

    return assembly


formpath = '/Users/mricardo/compas_dev/me/pattern/singular/crossvault/mesh-A2.json'

# formpath = '/Users/mricardo/compas_dev/compas_tno/data/input-form.json'
assemblypath = '/Users/mricardo/compas_dev/compas_tno/data/assembly-tno.json'
thk = 0.50

form = Mesh.from_json(formpath)

# Make dual
dual = extended_dual(form, Mesh)

# Make assembly
assembly_all = offset_dual(dual, thk)

# Remove last assembly block (Not sure why it is duplicated)
assembly = Assembly()
lastblock = len(list(assembly_all.blocks())) - 1
for node in assembly_all.nodes():
    block = assembly_all.node_block(node)
    if node < lastblock:
        assembly.add_block(block)

supports = [91, 101, 111, 81]

view = App()
view.add(form, show_faces=False, color=Color.red())
for node in assembly.nodes():
    block = assembly.node_block(node)
    if node in supports:
        view.add(block, opacity=0.5, color=Color.red())
    else:
        view.add(block, opacity=0.5)
view.show()

# assembly.to_json(assemblypath)
