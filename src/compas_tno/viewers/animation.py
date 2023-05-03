import time
import json

from compas.geometry import centroid_points
from compas.geometry import subtract_vectors

from compas_tno.plotters import TNOPlotter


def animation_from_optimisation(analysis, show_force=False, settings=None, record=False, interval=100, jump_each=1):
    """Make a 3D animated plot with the optimisation steps.

    Parameters
    ----------
    form : FormDiagram
        Form Diagram in which the animation should be based
    file_Xform : str
        Path for the ``.json`` file with the notal position of the form diagram during iterations
    force : ForceDiagram, optional
        Force diagram in whch the animation should be based, by default None
    file_Xforce : str, optional
        Path for the ``.json`` file with the notal position of the form diagram during iterations, by default None
    shape : Shape, optional
        The shape to plot in the animation, by default None
    settings : dict, optional
        Dictionary with settings to modify the viewer, by default None
    record : bool, optional
        Whether or not the animation should be saved as ``.gif``, by default False
    interval : int, optional
        Time interval (in ms) among iteration frames, by default 100
    jump_each : int, optional
        Interval tot jump among iteration frames, by default 1, in which all frames are shown
    """

    from compas_tno.algorithms import reciprocal_from_form
    from compas_tno.viewers import Viewer
    import compas_tno

    form = analysis.form
    # optimiser = analysis.optimiser
    shape = analysis.shape
    force = None
    if show_force:
        force = reciprocal_from_form(form)
    file_Xform = compas_tno.get('Xform.json')  # analysis.optimiser.Xform
    file_Xforce = compas_tno.get('Xforce.json')  # analysis.optimiser.Xforce

    viewer = Viewer(form, shape=shape)
    if settings:
        viewer.settings = settings
    viewer.initiate_app()

    with open(file_Xform, mode='r', encoding='utf-8') as f:
        Xform = json.load(f)

    if force:
        with open(file_Xforce, mode='r', encoding='utf-8') as f:
            Xforce = json.load(f)

    iterations_total = len(Xform)
    iterations = round(iterations_total / jump_each)

    print('Displaying {0} iterations from a total of {1} iterations'.format(iterations, iterations_total))

    out = None
    if record:
        if isinstance(record, str):
            out = record
        else:
            out = compas_tno.get('out.gif')
        print('Save gif at:', out)

    @viewer.app.on(interval=interval, frames=iterations, record=record, record_path=out)
    def update(f):

        print(f)

        if f == 1:
            time.sleep(5)

        viewer.clear()

        Xf = Xform[str(f * jump_each)]
        index = 0
        for vertex in form.vertices():
            viewer.thrust.vertex_attribute(vertex, 'x', Xf[index][0])
            viewer.thrust.vertex_attribute(vertex, 'y', Xf[index][1])
            viewer.thrust.vertex_attribute(vertex, 'z', Xf[index][2])
            index += 1

        viewer.draw_thrust()
        # viewer.draw_mesh(viewer.thrust)
        # viewer.draw_cracks()
        # viewer.draw_reactions()
        # viewer.draw_loads()
        # viewer.draw_shape()
        # viewer.draw_reactions()

        if show_force:
            _Xf = Xforce[str(f * jump_each)]
            index = 0
            for vertex in force.vertices():
                force.vertex_attribute(vertex, 'x', _Xf[index][0])
                force.vertex_attribute(vertex, 'y', _Xf[index][1])
                force.vertex_attribute(vertex, 'z', 0.0)
                index += 1
            force_centroid = centroid_points(force.vertices_attributes('xyz'))
            force_centroid0 = [-6.991104169236619, 66.12107435370979, 0.0]  # by hand because it isn't kept as a variable (it's deleted @ each passage)
            dc = subtract_vectors(force_centroid, force_centroid0)
            for vertex in force.vertices():
                x, y, _ = force.vertex_coordinates(vertex)
                force.vertex_attribute(vertex, 'x', x - dc[0])
                force.vertex_attribute(vertex, 'y', y - dc[1])
            viewer.draw_force(force)

    viewer.app.run()

    return


def animation_from_section(form, file_Xform):
    """Make a 3D animated plot with the xz section of the solution.

    Parameters
    ----------
    form : FormDiagram
        Form Diagram in which the animation should be based
    file_Xform : str
        Path for the ``.json`` file with the notal position of the form diagram during iterations

    """

    with open(file_Xform, mode='r', encoding='utf-8') as f:
        Xform = json.load(f)

    iterations = len(Xform)
    print('number of iterations', iterations)

    for i in range(iterations):

        print('Iteration:', i)

        Xi = Xform[str(i)]
        index = 0
        for vertex in form.vertices():
            form.vertex_attribute(vertex, 'x', Xi[index][0])
            form.vertex_attribute(vertex, 'y', Xi[index][1])
            form.vertex_attribute(vertex, 'z', Xi[index][2])
            index += 1

        plotter = TNOPlotter()
        meshartist = plotter.add(form)
        meshartist.draw_edges()
        plotter.show()

    return
