from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MultiMeshViewer
from compas.datastructures import Mesh

def view_intrados(shape):

    mesh = shape.intrados
    viewer = MeshViewer()
    viewer.mesh = mesh

    return viewer

def view_extrados(shape):

    mesh = shape.extrados
    viewer = MeshViewer()
    viewer.mesh = mesh


    return viewer

def view_middle(shape):

    mesh = shape.middle
    viewer = MeshViewer()
    viewer.mesh = mesh

    return viewer

def view_shapes(shape, cut_negatives=True):

    middle = shape.middle
    intrados = shape.intrados
    extrados = shape.extrados

    # Update to Delete the negatives...

    viewer = MultiMeshViewer()
    viewer.meshes = [middle, intrados, extrados]

    return viewer
