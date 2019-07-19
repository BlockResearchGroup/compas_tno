from compas_plotters import MeshPlotter

from numpy import array

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

def plot_form(form, radius=0.1, fix_width=False, max_width=10, simple=False):

    """ Extended load-path plotting of a FormDiagram

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

    q = [attr['q'] for u, v, attr in form.edges(True)]
    qmax  = max(abs(array(q)))
    lines = []

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), 'q')
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
    # round(form.get_vertex_attribute(i, 'pz'), 2)
    if radius:
        plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
        radius=radius, text={i: i for i in form.vertices()})
    plotter.draw_lines(lines)

    return plotter

def plot_force(force, radius=0.1, fix_width=False, max_width=10, simple=False):

    """ Extended load-path plotting of a ForceDiagram

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
    for key in force.edges():
        u,v = key
        lengths[key] = '{0:.2f}'.format(force.get_form_edge_attribute(form,key,'q'))

    plotter.draw_vertices(facecolor=vertexcolor, radius=radius)

    color = {key: '#00ff00' for key in force.edges() if force.get_form_edge_attribute(form, key, 'is_external')}
    width = {key: 2.0 for key in force.edges() if force.get_form_edge_attribute(form, key, 'is_external')}

    plotter.draw_edges(color=color, width=width, text=lengths)

    return plotter

def plot_grad(form, radius=0.1, fix_width=False, max_width=10, simple=True):

    """ Extended load-path plotting of a FormDiagram

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

