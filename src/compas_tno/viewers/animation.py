from compas.datastructures import Mesh
from compas_view2 import app
import time
import json


__all__ = [
    'animation_from_optimisation'
]


def animation_from_optimisation(form, file_Xform, force=None, file_Xforce=None, formscaling=None, forcescaling=None, shape=False, densities=False):
    """Make a 3D animated plot with the optimisation steps.
    """

    viewer = app.App()

    v, f = form.to_vertices_and_faces()
    form = Mesh.from_vertices_and_faces(v, f)
    key_index = form.key_index()
    obj = viewer.add(form)

    with open(file_Xform, mode='r', encoding='utf-8') as f:
        Xform = json.load(f)

    if force:
        v, f = force.to_vertices_and_faces()
        force = Mesh.from_vertices_and_faces(v, f)
        _key_index = force.key_index()
        _obj = viewer.add(force)

        with open(file_Xforce, mode='r', encoding='utf-8') as f:
            Xforce = json.load(f)

    if formscaling:
        print('No form scale available')
    if forcescaling:
        print('No force scale available')

    iterations = len(Xform)

    @viewer.on(interval=100, frames=iterations - 1)
    def update(f):
        print(f)
        if f == 0:
            time.sleep(5)

        Xf = Xform[str(f)]
        for vertex in form.vertices():
            index = key_index[vertex]
            form.vertex_attribute(vertex, 'x', Xf[index][0])
            form.vertex_attribute(vertex, 'y', Xf[index][1])
            form.vertex_attribute(vertex, 'z', Xf[index][2])
        # Do transformation if formscale
        obj.update()

        if force:
            _Xf = Xforce[str(f)]
            for vertex in force.vertices():
                index = _key_index[vertex]
                force.vertex_attribute(vertex, 'x', _Xf[index][0])
                force.vertex_attribute(vertex, 'y', _Xf[index][1])
                force.vertex_attribute(vertex, 'z', _Xf[index][2])
        # Do transformation if forcescale
            _obj.update()

    viewer.run()

    return
