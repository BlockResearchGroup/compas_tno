from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas_tna.diagrams import FormDiagram

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'plot_force',
    'plot_dual',
]


def plot_force(force, form, show_length=False, radius=0.1, fix_width=False, max_width=10, simple=False, color_inds=True):
    """ Extended plotting of a ForceDiagram

    Parameters
    ----------
    force : obj
        ForceDiagram to plot.
    form : obj
        FormDiagram realated to plot.
    show_length : bool
        Show the length of edges in the force diagram, that relate to the magnitude of the force in the FormDiagram.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : float
        Maximum edge width.
    simple : bool
        Simple red and blue colour plotting.
    color_inds : bool
        Color the independent edges.

    Returns
    -------
    obj
        Plotter object.

    """

    plotter = MeshPlotter(force, figsize=(12, 8), tight=True)

    vertexcolor = {key: (1.0, 0.9, 0.9) for key in force.vertices() if not form.face_attribute(key, '_is_loaded')}

    radius = {key: 0.05 for key in force.vertices()}
    radius.update({key: 0.1 for key in force.vertices() if not form.face_attribute(key, '_is_loaded')})

    lengths = {}
    color = {}
    width = {}
    for key in force.edges():
        lengths[key] = '{0:.2f}'.format(force.form_edge_attribute(form, key, 'q'))
        if force.form_edge_attribute(form, key, '_is_external'):
            color[key] = '#00ff00'
            width[key] = 2.0
        if color_inds:
            if force.form_edge_attribute(form, key, 'is_ind'):
                color[key] = '#ff0000'
                width[key] = 3.0
            pass

    plotter.draw_vertices(facecolor=vertexcolor, radius=radius)
    if show_length:
        plotter.draw_edges(color=color, width=width, text=lengths)
    else:
        plotter.draw_edges(color=color, width=width)

    return plotter


def plot_dual(form):
    """ Plot the dual of a FormDiagram

    Parameters
    ----------
    form : obj
        FormDiagram to plot the dual from.

    Returns
    -------
    obj
        Plotter object.

    """

    lines = []
    dual = FormDiagram.dual(form, Mesh)

    for u, v in dual.edges():
        lines.append({
            'start': dual.vertex_coordinates(u, 'xy'),
            'end': dual.vertex_coordinates(v, 'xy'),
            'color': '#000000',
            'width': 0.5
        })
    for u, v in form.edges():
        lines.append({
            'start': form.vertex_coordinates(u, 'xy'),
            'end': form.vertex_coordinates(v, 'xy'),
            'color': '#FF0000',
            'width': 1.0
        })

    plotter = MeshPlotter(dual, figsize=(10, 6))
    plotter.draw_lines(lines)

    return plotter
