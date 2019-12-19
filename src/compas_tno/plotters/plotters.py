from compas_plotters import MeshPlotter

from numpy import array

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from math import sqrt

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'plot_form',
    'plot_force',
    'plot_grad',
    'plot_dual',
    'plot_form_xz',
    'plot_form_joints',
]


def plot_form(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=True, thick='q', heights=False, show_edgeuv=False, save=None):
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
    ----------
    obj
        Plotter object.

    """

    uv_i = form.uv_index()
    q = [attr[thick] for u, v, attr in form.edges(True)]
    qmax = max(abs(array(q)))
    lines = []
    i = 0

    for u, v in form.edges_where({'is_edge': True}):
        qi = form.get_edge_attribute((u, v), thick)
        l = form.edge_length(u, v)
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
            # text = str(u) + ',' + str(v)
            text = str(i)
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

        i = i + 1

    rad_colors = {}
    for key in form.vertices_where({'is_fixed': True}):
        rad_colors[key] = '#aaaaaa'
    for key in form.vertices_where({'rol_x': True}):
        rad_colors[key] = '#ffb733'
    for key in form.vertices_where({'rol_y': True}):
        rad_colors[key] = '#ffb733'

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        if heights:
            plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
                                  radius=radius, text={i: i for i in form.vertices()})  # form.get_vertex_attribute(i, 'z')
        else:
            plotter.draw_vertices(facecolor=rad_colors, radius=radius)

    plotter.draw_lines(lines)
    if save:
        plotter.save(save)

    return plotter


def plot_force(force, form, show_length=False, radius=0.1, fix_width=False, max_width=10, simple=False, color_inds=True):
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
        u, v = key
        lengths[key] = '{0:.2f}'.format(force.get_form_edge_attribute(form, key, 'q'))
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
    qmax = max(abs(array(q)))
    lines = []

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), 'dq')
        l = form.edge_length(u, v)
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
            'text': round(qi, 2),
        })

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})}, radius=radius)
    plotter.draw_lines(lines)

    return plotter


def plot_dual(form):

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


def plot_form_xz(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=True, thick='q', heights=False, show_edgeuv=False, save=None, thk=0.20, plot_reactions=False, joints=False, cracks=False):
    """ Plor of a 2D diagrma in the XZ plane

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

    i_k = form.index_key()
    q = [attr[thick] for u, v, attr in form.edges(True)]
    qmax = max(abs(array(q)))
    lines = []
    xs = []
    reac_lines = []

    for key in form.vertices():
        xs.append(form.vertex_coordinates(key)[0])
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            x, _, z = form.vertex_coordinates(key)
            if z > 0.0:
                rz = abs(form.get_vertex_attribute(key, 'rz'))
                rx = form.get_vertex_attribute(key, 'rx')
                reac_line = [x, z, x + z * rx / rz, 0.0]
                reac_lines.append(reac_line)
                # reac_x.append(x + z*rx/rz)
                # reac_z.append(0.0)
                # reac_x.append(x)
                # reac_z.append(z)
                # rx attributes not right!!

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), thick)
        l = form.edge_length(u, v)
        uv_i = form.uv_index

        if simple:
            if qi > 0:
                colour = ['00', '00', '00']
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
            'start': [form.vertex_coordinates(u)[0], form.vertex_coordinates(u)[2]],
            'end':   [form.vertex_coordinates(v)[0], form.vertex_coordinates(v)[2]],
            'color': 'FF0000',
            'width': width,
            'text': text,
        })

    try:
        Re = form.attributes['Re']
        Ri = form.attributes['Ri']
    except:
        Re = 1.20  # (max(xs) - min(xs))/2 + thk/2
        Ri = 1.00  # (max(xs) - min(xs))/2 - thk/2

    xc = sum(xs)/len(xs)
    discr = 200
    # print('Visualisation on Re: {0:.3f} / Ri: {1:.3f}'.format(Re,Ri))

    for R in [Re, Ri]:
        for i in range(discr):
            lines.append({
                'start': [xc-R+2*R*i/discr, sqrt(abs(R**2 - (2*R*i/discr-R)**2))],
                'end':   [xc-R+2*R*(i+1)/discr, sqrt(abs(R**2 - (2*R*(i+1)/discr-R)**2))],
                'color': '000000',
                'width': 0.5,
            })

    if plot_reactions:
        for reac_line in reac_lines:
            lines.append({
                'start': [reac_line[0], reac_line[1]],
                'end':   [reac_line[2], reac_line[3]],
                'color': ''.join(colour),
                'width': width,
            })

    if joints:
        joints = form.attributes['joints']
        for i in joints:
            lines.append({
                'start': [joints[i][0][0], joints[i][0][2]],
                'end':   [joints[i][1][0], joints[i][1][2]],
                'color': '000000',
                'width': 0.25,
            })

    vertices = []
    if cracks:
        cracks_lb, cracks_ub = form.attributes['cracks']
        for i in cracks_ub:
            key = i_k[i]
            x, _, _ = form.vertex_coordinates(key)
            z = form.get_vertex_attribute(key, 'ub')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })
        for i in cracks_lb:
            key = i_k[i]
            x, _, _ = form.vertex_coordinates(key)
            z = form.get_vertex_attribute(key, 'lb')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })

    nodes = []
    if radius:
        for key in form.vertices():
            x, _, z = form.vertex_coordinates(key)
            if form.get_vertex_attribute(key, 'is_fixed') is True:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '000000',
                    'facecolor': 'aaaaaa',
                })
            if abs(form.get_vertex_attribute(key, 'ub') - z) < 1e-5:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': 'FF0000',
                    'facecolor': 'FF0000',
                })
            if abs(form.get_vertex_attribute(key, 'lb') - z) < 1e-5:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '0000FF',
                    'facecolor': '0000FF',
                })

    plotter = MeshPlotter(form, figsize=(10, 10))
    # round(form.get_vertex_attribute(i, 'pz'), 2)
    # if radius:
    #     if heights:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius, text={i: i for i in form.vertices()}) # form.get_vertex_attribute(i, 'z')
    #     else:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius)

    # plotter.draw_vertices(radius= {i : form.get_vertex_attribute(i, 'px')/100 for i in form.vertices()}) # form.get_vertex_attribute(i, 'z')

    plotter.draw_lines(lines)
    plotter.draw_points(vertices)
    plotter.draw_points(nodes)

    if save:
        plotter.save(save)

    return plotter


def plot_form_joints(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=True, thick='q', heights=False, show_edgeuv=False, save=None, thk=0.20, plot_reactions=False):
    """ Plor of a 2D diagrma in the XZ plane

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
    qmax = max(abs(array(q)))
    lines = []
    xs = []
    reac_lines = []

    for key in form.vertices():
        xs.append(form.vertex_coordinates(key)[0])
        if form.get_vertex_attribute(key, 'is_fixed') == True:
            x, _, z = form.vertex_coordinates(key)
            if z > 0.0:
                rz = abs(form.get_vertex_attribute(key, 'rz'))
                rx = form.get_vertex_attribute(key, 'rx')
                reac_line = [x, z, x + z*rx/rz, 0.0]
                reac_lines.append(reac_line)
                # reac_x.append(x + z*rx/rz)
                # reac_z.append(0.0)
                # reac_x.append(x)
                # reac_z.append(z)

    for u, v in form.edges():
        qi = form.get_edge_attribute((u, v), thick)
        l = form.edge_length(u, v)
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
            'start': [form.vertex_coordinates(u)[0], form.vertex_coordinates(u)[2]],
            'end':   [form.vertex_coordinates(v)[0], form.vertex_coordinates(v)[2]],
            'color': 'FF0000',
            'width': width,
            'text': text,
        })

    try:
        Re = form.attributes['Re']
        Ri = form.attributes['Ri']
    except:
        Re = 1.20  # (max(xs) - min(xs))/2 + thk/2
        Ri = 1.00  # (max(xs) - min(xs))/2 - thk/2

    xc = sum(xs)/len(xs)
    discr = 200
    # print('Visualisation on Re: {0:.3f} / Ri: {1:.3f}'.format(Re,Ri))

    for R in [Re, Ri]:
        for i in range(discr):
            lines.append({
                'start': [xc-R+2*R*i/discr, sqrt(abs(R**2 - (2*R*i/discr-R)**2))],
                'end':   [xc-R+2*R*(i+1)/discr, sqrt(abs(R**2 - (2*R*(i+1)/discr-R)**2))],
                'color': '000000',
                'width': 0.5,
            })

    if plot_reactions:
        for reac_line in reac_lines:
            lines.append({
                'start': [reac_line[0], reac_line[1]],
                'end':   [reac_line[2], reac_line[3]],
                'color': ''.join(colour),
                'width': width,
            })

    plotter = MeshPlotter(form, figsize=(10, 10))
    # round(form.get_vertex_attribute(i, 'pz'), 2)
    # if radius:
    #     if heights:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius, text={i: i for i in form.vertices()}) # form.get_vertex_attribute(i, 'z')
    #     else:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius)

    # plotter.draw_vertices(radius= {i : form.get_vertex_attribute(i, 'px')/100 for i in form.vertices()}) # form.get_vertex_attribute(i, 'z')

    plotter.draw_lines(lines)
    if save:
        plotter.save(save)

    return plotter
