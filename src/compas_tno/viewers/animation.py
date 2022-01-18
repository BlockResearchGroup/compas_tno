from compas.datastructures import Mesh
from compas_view2 import app
import time
import json

from compas_tno.viewers.viewer import Viewer


__all__ = [
    'animation_from_optimisation',
    'save_geometry_at_iterations'
]


def save_geometry_at_iterations(form, optimiser, shape=None, force=None):

    import json
    import compas_tno
    from numpy import array
    from compas_tno.algorithms import xyz_from_xopt
    from compas_tno.algorithms import reciprocal_from_form

    M = optimiser.M  # matrices of the problem

    file_qs = compas_tno.get('output.json')
    file_Xform = compas_tno.get('Xform.json')

    if force:
        file_Xforce = compas_tno.get('Xforce.json')
        Xforce = {}

    with open(file_qs, mode='r', encoding='utf-8') as f:
        data = json.load(f)

    Xform = {}

    iterations = len(data['iterations'])

    for i in range(iterations):
        xopt_i = array(data['iterations'][str(i)]).reshape(-1, 1)
        M = xyz_from_xopt(xopt_i, M)
        Xform_i = M.X.tolist()
        Xform[str(i)] = Xform_i

        if force:
            j = 0
            for key in form.vertices():
                form.vertex_attributes(key, 'xyz', M.X[j].tolist())
                j += 1
            k = 0
            for edge in form.edges_where({'_is_edge': True}):
                form.edge_attribute(edge, 'q', M.q.flatten()[k])
                k += 1

            force = reciprocal_from_form(form)
            Xforce_i = force.vertices_attributes('xyz')
            Xforce[str(i)] = Xforce_i

    with open(file_Xform, mode='w', encoding='utf-8') as f:
        json.dump(Xform, f)

    print('Form Geometry saved @:', file_Xform)

    if force:
        with open(file_Xforce, mode='w', encoding='utf-8') as f:
            json.dump(Xforce, f)

        print('Force Geometry saved @:', file_Xforce)

    return


def animation_from_optimisation(form, file_Xform, force=None, file_Xforce=None, forcescaling=None, forcetranslation=None, shape=False, densities=False, interval=100):
    """ Make a 3D animated plot with the optimisation steps.
    """

    viewer = Viewer(form)

    v, f = form.to_vertices_and_faces()
    form = Mesh.from_vertices_and_faces(v, f)
    key_index = form.key_index()

    with open(file_Xform, mode='r', encoding='utf-8') as f:
        Xform = json.load(f)

    if force:
        with open(file_Xforce, mode='r', encoding='utf-8') as f:
            Xforce = json.load(f)

    if forcescaling:
        viewer.settings['force.scale'] = forcescaling
    if forcetranslation:
        viewer.settings['force.translation'] = forcetranslation

    iterations = len(Xform)
    print('number of iterations', iterations)

    @viewer.app.on(interval=interval, frames=iterations)
    def update(f):

        print(f)

        if f == 0:
            time.sleep(5)

        viewer.clear()

        Xf = Xform[str(f)]
        for vertex in form.vertices():
            index = key_index[vertex]
            viewer.thrust.vertex_attribute(vertex, 'x', Xf[index][0])
            viewer.thrust.vertex_attribute(vertex, 'y', Xf[index][1])
            viewer.thrust.vertex_attribute(vertex, 'z', Xf[index][2])

        viewer.view_thrust()
        viewer.view_cracks()
        viewer.view_shape()
        viewer.view_reactions()

        if force:
            _Xf = Xforce[str(f)]
            index = 0
            for vertex in force.vertices():
                force.vertex_attribute(vertex, 'x', _Xf[index][0])
                force.vertex_attribute(vertex, 'y', _Xf[index][1])
                force.vertex_attribute(vertex, 'z', _Xf[index][2])
                index += 1
            viewer.view_force(force)

    viewer.app.run()

    return
