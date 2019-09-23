from compas_plotters import MeshPlotter

from numpy import array

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'plot_form',
    'plot_force',
    'plot_grad',
    'plot_dual',
]

def plot_form(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q =True, thick = 'q', heights = False, show_edgeuv=False):

    """ Extended plotting of a FormDiagram

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : float
        Maximum edge width.
    simple : bool
        Simple red and blue colour plotting.

    Returns
    -------
    obj
        Plotter object.

    """

    q = [attr[thick] for u, v, attr in form.edges(True)]
    qmax  = 1 # max(abs(array(q))) # 
    lines = []

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), thick)
        l = form.edge_length(u,v)
        uv_i = form.uv_index

        if simple:
            if qi > 0:
                colour = ['ff', '00', '00']
            elif qi < 0:
                colour = ['00', '00', 'ff']
            else:
                colour = ['aa', 'aa', 'aa']

        else:
            colour = ['ff', '00', '00']
            if qi > 0:
                colour[0] = 'ff'
            if form.get_edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.get_edge_attribute((u, v), 'is_ind'):
                # colour[2] = 'ff'
                colour[0] = '00'
                colour[2] = '80'

        width = max_width if fix_width else (qi / qmax) * max_width

        
        if show_edgeuv:
            text = str(u) + ',' + str(v)
        elif show_q:
            text = round(qi, 2)
        else:
            text = ''

        lines.append({
            'start': form.vertex_coordinates(u),
            'end':   form.vertex_coordinates(v),
            'color': ''.join(colour),
            'width': width,
            'text': text,
        })

    plotter = MeshPlotter(form, figsize=(10, 10))
    # round(form.get_vertex_attribute(i, 'pz'), 2)
    if radius:
        if heights:
            plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
            radius=radius, text={i: form.get_vertex_attribute(i, 'pz') for i in form.vertices()})
        else:
            plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
            radius=radius)

    plotter.draw_lines(lines)

    return plotter

def plot_force(force, form, show_length = False, radius=0.1, fix_width=False, max_width=10, simple=False, color_inds=True):

    """ Extended plotting of a Formdiagram

    Parameters
    ----------
    form : obj
        ForceDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : float
        Maximum edge width.
    simple : bool
        Simple red and blue colour plotting.

    Returns
    -------
    obj
        Plotter object.

    """


    plotter = MeshPlotter(force, figsize=(12, 8), tight=True)

    vertexcolor = {key: (1.0, 0.9, 0.9) for key in force.vertices() if not form.get_face_attribute(key, 'is_loaded')}

    radius = {key: 0.05 for key in force.vertices()}
    radius.update({key: 0.1 for key in force.vertices() if not form.get_face_attribute(key, 'is_loaded')})

    lengths = {}
    color = {}
    width = {}
    for key in force.edges():
        u,v = key
        lengths[key] = '{0:.2f}'.format(force.get_form_edge_attribute(form,key,'q'))
        if force.get_form_edge_attribute(form, key, 'is_external'):
            color[key] = '#00ff00' 
            width[key] = 2.0
        if color_inds:
            if force.get_form_edge_attribute(form, key, 'is_ind'):
                color[key] = '#ff0000' 
                width[key] = 3.0
            pass


    plotter.draw_vertices(facecolor=vertexcolor, radius=radius)
    if show_length:
        plotter.draw_edges(color=color, width=width, text=lengths)
    else:
        plotter.draw_edges(color=color, width=width)

    return plotter

def plot_grad(form, radius=0.1, fix_width=False, max_width=10, simple=True):

    """ Extended load-path plotting of a the Gradient

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : float
        Maximum edge width.
    simple : bool
        Simple red and blue colour plotting.

    Returns
    -------
    obj
        Plotter object.

    """

    q = [attr['dq'] for u, v, attr in form.edges(True)]
    qmax  = max(abs(array(q)))
    lines = []

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), 'dq')
        l = form.edge_length(u,v)
        uv_i = form.uv_index

        if simple:
            if qi > 0:
                colour = ['ff', '00', '00']
            elif qi < 0:
                colour = ['00', '00', 'ff']
            else:
                colour = ['aa', 'aa', 'aa']

        else:
            colour = ['00', '00', '00']
            if qi > 0:
                colour[0] = 'ff'
            if form.get_edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.get_edge_attribute((u, v), 'is_ind'):
                colour[2] = 'ff'

        width = max_width if fix_width else (qi / qmax) * max_width

        lines.append({
            'start': form.vertex_coordinates(u),
            'end':   form.vertex_coordinates(v),
            'color': ''.join(colour),
            'width': width,
            'text':round(qi, 2),
        })

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})}, radius=radius)
    plotter.draw_lines(lines)

    return plotter


def plot_dual(form):

    lines = []
    dual = FormDiagram.dual(form,Mesh)

    for u, v in dual.edges():
        lines.append({
            'start': dual.vertex_coordinates(u, 'xy'),
            'end'  : dual.vertex_coordinates(v, 'xy'),
            'color': '#000000',
            'width': 0.5
        })
    for u, v in form.edges():
        lines.append({
            'start': form.vertex_coordinates(u, 'xy'),
            'end'  : form.vertex_coordinates(v, 'xy'),
            'color': '#FF0000',
            'width': 1.0
        })

    plotter = MeshPlotter(dual, figsize=(10, 6))
    plotter.draw_lines(lines)

    return plotter
