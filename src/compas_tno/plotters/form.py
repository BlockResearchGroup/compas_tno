import compas_tno
from compas_plotters import MeshPlotter
from numpy import array
from numpy import linspace

import matplotlib.pyplot as plt
from matplotlib import cm

from compas_tna.diagrams import FormDiagram
from compas.utilities import geometric_key
from compas_plotters import Plotter
from math import sqrt
import math
import os


__all__ = [
    'plot_form',
    'plot_form_xz',
    'plot_forms_xz',
    'plot_independents',
    'plot_gif_forms_xz',
    'plot_sym_inds',
]


def plot_form(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=True, thick='q', heights=False, show_edgeuv=False, cracks=False, save=None):
    """ Extended plotting of a FormDiagram

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float (0.05)
        Radius of vertex markers.
    fix_width : bool (False)
        Fix edge widths as constant.
    max_width : bool (10)
        Maximum width of the plot.
    max_width : float (False)
        Maximum edge width.
    simple : bool (True)
        Simple red and blue colour plotting.
    show_q : bool (True)
        Show the force densities on the edges.
    thick : str ('q')
        Attribute that the thickness of the form should be related to.
    heights : bool (False)
        Plot the heights of the nodes.
    show_edgeuv : bool (False)
        Show u,v of the edges.
    cracks : bool (False)
        If true highlight the location of the nodes touching intrados (blue) and extrados (green).
    save : str (None)
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """

    uv_i = form.uv_index()
    q = [form.edge_attribute((u,v), thick) for u, v in form.edges_where({'_is_edge': True})]
    qmax = max(abs(array(q)))
    lines = []
    i = 0

    for u, v in form.edges_where({'_is_edge': True}):
        qi = form.edge_attribute((u, v), thick)
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
            if form.edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.edge_attribute((u, v), 'is_ind'):
                # colour[2] = 'ff'
                colour[0] = '00'
                colour[2] = '80'

        width = max_width if fix_width else (qi / qmax) * max_width

        if show_edgeuv:
            text = str(u) + ',' + str(v)
            # text = str(i)
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
    if cracks:
        for key in form.vertices():
            ub = form.vertex_attribute(key, 'ub')
            lb = form.vertex_attribute(key, 'lb')
            z = form.vertex_attribute(key, 'z')
            if abs(z - ub) < 10e-4:
                rad_colors[key] = '#008000' # Green extrados
            elif abs(z - lb) < 10e-4:
                rad_colors[key] = '#0000FF' # Blue intrados
            elif z - ub > 0 or lb - z > 0:
                rad_colors[key] = '#000000' # Black outside



    for key in form.vertices_where({'is_fixed': True}):
        rad_colors[key] = '#aaaaaa'
    for key in form.vertices_where({'rol_x': True}):
        rad_colors[key] = '#ffb733'
    for key in form.vertices_where({'rol_y': True}):
        rad_colors[key] = '#ffb733'

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        plotter.draw_vertices(facecolor=rad_colors, radius=radius)
        if heights:
            plotter.draw_vertices(keys= [i for i in form.vertices_where({'is_fixed': True})], facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
                                  radius=radius, text={i: [round(form.vertex_attribute(i, 'lb'),3),round(form.vertex_attribute(i, 'ub'),3),round(form.vertex_attribute(i, 'z'),3)] for i in form.vertices()})  # form.vertex_attribute(i, 'z')

    plotter.draw_lines(lines)
    if save:
        plotter.save(save)

    return plotter


def plot_form_xz(form, shape, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=False, plot_reactions=True, cracks=False, save=False):
    """ Plot a FormDiagram in axis xz

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    shape: obj
        Shape to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : bool
        Maximum width of the plot.
    max_width : float
        Maximum edge width.
    simple : bool
        Simple red and blue colour plotting.
    show_q : bool
        Show the force densities on the edges.
    thick : str
        Attribute that the thickness of the form should be related to.
    heights : bool
        Plot the heights of the nodes.
    show_edgeuv : bool
        Show u,v of the edges.
    save : str
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """
    i_k = form.index_key()
    q = [form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]
    qmax = max(abs(array(q)))
    lines = []
    xs = []
    reac_lines = []
    discr = 100

    for key in form.vertices():
        xs.append(form.vertex_coordinates(key)[0])
        if form.vertex_attribute(key, 'is_fixed') == True:
            x, _, z = form.vertex_coordinates(key)
            if z > 0.0:
                rz = abs(form.vertex_attribute(key, '_rz'))
                rx = form.vertex_attribute(key, '_rx')
                reac_line = [x, z, x + z * rx / rz, 0.0]
                reac_lines.append(reac_line)


    for u, v in form.edges():
        qi = form.edge_attribute((u, v), 'q')

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
            if form.edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.edge_attribute((u, v), 'is_ind'):
                # colour[2] = 'ff'
                colour[0] = '00'
                colour[2] = '80'

        width = max_width if fix_width else (qi / qmax) * max_width

        if show_q:
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

    if shape.data['type'] == 'arch':

        H = shape.data['H']
        L = shape.data['L']
        thk = shape.data['thk']
        R = H / 2 + (L**2 / (8 * H))
        zc = R - H
        re = R + thk/2
        ri = R - thk/2
        x = sqrt(re**2 - zc**2)
        spr_e = math.acos(zc/re)
        # tot_angle_i = 2*spr_i
        tot_angle_e = 2*spr_e
        # angle_init_i = (math.pi - tot_angle_i)/2
        angle_init_e = (math.pi - tot_angle_e)/2
        # an_i = tot_angle_i / discr
        an_e = tot_angle_e / discr
        xc = L/2

        for i in range(discr):
            angle_i = angle_init_e + i * an_e
            angle_f = angle_init_e + (i + 1) * an_e
            for r_ in [ri, re]:
                xi = xc - r_ * math.cos(angle_i)
                xf = xc - r_ * math.cos(angle_f)
                zi = r_ * math.sin(angle_i) - zc
                zf = r_ * math.sin(angle_f) - zc
                lines.append({
                    'start': [xi, zi],
                    'end':   [xf, zf],
                    'color': '000000',
                    'width': 0.5,
                })

        if plot_reactions:
            for reac_line in reac_lines:
                lines.append({
                    'start': [reac_line[0], reac_line[1]],
                    'end':   [reac_line[2], reac_line[3]],
                    'color': ''.join(['00', '00', '00']),
                    'width': max_width,
                })

        vertices = []
        if cracks:
            cracks_lb, cracks_ub = form.attributes['cracks']
            for i in cracks_ub:
                key = i_k[i]
                x, _, _ = form.vertex_coordinates(key)
                z = form.vertex_attribute(key, 'ub')
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'color': '000000',
                })
            for i in cracks_lb:
                key = i_k[i]
                x, _, _ = form.vertex_coordinates(key)
                z = form.vertex_attribute(key, 'lb')
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'color': '000000',
                })

        nodes = []
        if radius:
            for key in form.vertices():
                x, _, z = form.vertex_coordinates(key)
                if form.vertex_attribute(key, 'is_fixed') is True:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '000000',
                        'facecolor': 'aaaaaa',
                    })
                if abs(form.vertex_attribute(key, 'ub') - z) < 1e-3:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '008000',
                        'facecolor': '008000',
                    })
                if abs(form.vertex_attribute(key, 'lb') - z) < 1e-3:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '0000FF',
                        'facecolor': '0000FF',
                    })

    plotter = MeshPlotter(form, figsize=(10, 10))
    plotter.draw_lines(lines)
    plotter.draw_points(vertices)
    plotter.draw_points(nodes)

    if save:
        plotter.save(save)

    return plotter


def plot_forms_xz(forms, shape, radius=0.05, fix_width=False, max_width=10, plot_reactions=True, cracks=False, save=False):
    """ Plot multiple FormDiagrams in axis xz

    Parameters
    ----------
    form : list
        List of FormDiagrams to plot.
    shape: obj
        Shape to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : bool
        Maximum width of the plot.
    plot_reactions : bool
        Plot reactions.
    cracks : bool
        Show cracks - points where the thrust line touches the masonry limits.
    save : str
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """

    plotter = Plotter()
    vertices = []
    lines = []
    xs = []
    reac_lines = []
    nodes = []

    if shape.data['type'] == 'arch':
        discr = 100
        H = shape.data['H']
        L = shape.data['L']
        thk = shape.data['thk']
        R = H / 2 + (L**2 / (8 * H))
        zc = R - H
        re = R + thk/2
        ri = R - thk/2
        x = sqrt(re**2 - zc**2)
        spr_e = math.acos(zc/re)
        tot_angle_e = 2*spr_e
        angle_init_e = (math.pi - tot_angle_e)/2
        an_e = tot_angle_e / discr
        xc = L/2

        for i in range(discr):
            angle_i = angle_init_e + i * an_e
            angle_f = angle_init_e + (i + 1) * an_e
            for r_ in [ri, re]:
                xi = xc - r_ * math.cos(angle_i)
                xf = xc - r_ * math.cos(angle_f)
                zi = r_ * math.sin(angle_i) - zc
                zf = r_ * math.sin(angle_f) - zc
                lines.append({
                    'start': [xi, zi],
                    'end':   [xf, zf],
                    'color': '000000',
                    'width': 0.5,
                })

    for form in forms:
        i_k = form.index_key()
        q = [form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]
        qmax = max(abs(array(q)))

        for key in form.vertices():
            xs.append(form.vertex_coordinates(key)[0])
            if form.vertex_attribute(key, 'is_fixed') == True:
                x, _, z = form.vertex_coordinates(key)
                if z > 0.0:
                    rz = abs(form.vertex_attribute(key, '_rz'))
                    rx = form.vertex_attribute(key, '_rx')
                    reac_line = [x, z, x + z * rx / rz, 0.0]
                    reac_lines.append(reac_line)

        for u, v in form.edges():
            qi = form.edge_attribute((u, v), 'q')
            width = max_width if fix_width else (qi / qmax) * max_width
            lines.append({
                'start': [form.vertex_coordinates(u)[0], form.vertex_coordinates(u)[2]],
                'end':   [form.vertex_coordinates(v)[0], form.vertex_coordinates(v)[2]],
                'color': 'FF0000',
                'width': width,
            })

        if plot_reactions:
            for reac_line in reac_lines:
                lines.append({
                    'start': [reac_line[0], reac_line[1]],
                    'end':   [reac_line[2], reac_line[3]],
                    'color': ''.join(['00', '00', '00']),
                    'width': max_width,
                })

        if cracks:
            cracks_lb, cracks_ub = form.attributes['cracks']
            for i in cracks_ub:
                key = i_k[i]
                x, _, _ = form.vertex_coordinates(key)
                z = form.vertex_attribute(key, 'ub')
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'color': '000000',
                })
            for i in cracks_lb:
                key = i_k[i]
                x, _, _ = form.vertex_coordinates(key)
                z = form.vertex_attribute(key, 'lb')
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'color': '000000',
                })
        if radius:
            for key in form.vertices():
                x, _, z = form.vertex_coordinates(key)
                if form.vertex_attribute(key, 'is_fixed') is True:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '000000',
                        'facecolor': 'aaaaaa',
                    })
                if abs(form.vertex_attribute(key, 'ub') - z) < 1e-3:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '008000',
                        'facecolor': '008000',
                    })
                if abs(form.vertex_attribute(key, 'lb') - z) < 1e-3:
                    nodes.append({
                        'pos': [x, z],
                        'radius': radius,
                        'edgecolor': '0000FF',
                        'facecolor': '0000FF',
                    })

    plotter.draw_lines(lines)
    plotter.draw_points(vertices)
    plotter.draw_points(nodes)

    if save:
        plotter.save(save)

    return plotter


def lines_and_points_from_form(form, plot_reactions, cracks, radius, max_width, fix_width, hide_negative=False):
    vertices = []
    lines = []
    xs = []
    reac_lines = []
    i_k = form.index_key()
    q = [form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]
    qmax = max(abs(array(q)))

    for key in form.vertices():
        xs.append(form.vertex_coordinates(key)[0])
        if form.vertex_attribute(key, 'is_fixed') == True:
            x, _, z = form.vertex_coordinates(key)
            if z > 0.0:
                rz = abs(form.vertex_attribute(key, '_rz'))
                rx = form.vertex_attribute(key, '_rx')
                reac_line = [x, z, x + z * rx / rz, 0.0]
                reac_lines.append(reac_line)

    for u, v in form.edges():
        qi = form.edge_attribute((u, v), 'q')
        width = max_width if fix_width else (qi / qmax) * max_width
        if not hide_negative or (form.vertex_coordinates(u)[2] > -10e-4 and form.vertex_coordinates(v)[2] > -10e-4):
            lines.append({
                'start': [form.vertex_coordinates(u)[0], form.vertex_coordinates(u)[2]],
                'end':   [form.vertex_coordinates(v)[0], form.vertex_coordinates(v)[2]],
                'color': 'FF0000',
                'width': width,
            })
        elif form.vertex_coordinates(u)[2] > -10e-4 or form.vertex_coordinates(v)[2]> -10e-4:
            m = (form.vertex_coordinates(u)[2] - form.vertex_coordinates(v)[2])/(form.vertex_coordinates(u)[0] - form.vertex_coordinates(v)[0])
            xzero = -1 * form.vertex_coordinates(u)[2]/m + form.vertex_coordinates(u)[0]
            b = max(form.vertex_coordinates(u)[2], form.vertex_coordinates(v)[2])
            a = form.vertex_coordinates(u)[0] if b == form.vertex_coordinates(u)[2] else form.vertex_coordinates(v)[0]
            lines.append({
                'start': [xzero, 0.0],
                'end':   [a, b],
                'color': 'FF0000',
                'width': width,
            })

    if plot_reactions:
        for reac_line in reac_lines:
            lines.append({
                'start': [reac_line[0], reac_line[1]],
                'end':   [reac_line[2], reac_line[3]],
                'color': ''.join(['00', '00', '00']),
                'width': max_width,
            })

    if cracks:
        cracks_lb, cracks_ub = form.attributes['cracks']
        for i in cracks_ub:
            key = i_k[i]
            x, _, _ = form.vertex_coordinates(key)
            z = form.vertex_attribute(key, 'ub')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })
        for i in cracks_lb:
            key = i_k[i]
            x, _, _ = form.vertex_coordinates(key)
            z = form.vertex_attribute(key, 'lb')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })
    if radius:
        for key in form.vertices():
            x, _, z = form.vertex_coordinates(key)
            if (form.vertex_attribute(key, 'is_fixed') is True and not hide_negative) or (hide_negative and form.vertex_attribute(key, 'is_fixed') is True and z > -10e-4):
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '000000',
                    'facecolor': 'aaaaaa',
                })
            if abs(form.vertex_attribute(key, 'ub') - z) < 1e-3:
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '008000',
                    'facecolor': '008000',
                })
            if abs(form.vertex_attribute(key, 'lb') - z) < 1e-3:
                vertices.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '0000FF',
                    'facecolor': '0000FF',
                })

    return lines, vertices


def plot_gif_forms_xz(forms, shape, radius=0.05, fix_width=False, max_width=10, plot_reactions=True, cracks=False, hide_negative=True, save=False, delay=0.5):
    """ Plot multiple FormDiagrams in axis xz in a gif.

    Parameters
    ----------
    form : list
        List of FormDiagrams to plot.
    shape: obj
        Shape to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : bool
        Maximum width of the plot.
    plot_reactions : bool
        Plot reactions.
    cracks : bool
        Show cracks - points where the thrust line touches the masonry limits.
    hide_negative : bool
        Hide negative points.
    save : str
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """
    plotter = Plotter(figsize=(10, 10))
    vertices = []
    lines_arch = []
    images = []
    tempfolder = compas_tno.get('/temp')
    img_count = 0
    pattern = 'image_{}.png'

    if shape.data['type'] == 'arch':
        discr = 100
        H = shape.data['H']
        L = shape.data['L']
        thk = shape.data['thk']
        R = H / 2 + (L**2 / (8 * H))
        zc = R - H
        re = R + thk/2
        ri = R - thk/2
        spr_e = math.acos(zc/re)
        tot_angle_e = 2*spr_e
        angle_init_e = (math.pi - tot_angle_e)/2
        an_e = tot_angle_e / discr
        xc = L/2

        for i in range(discr):
            angle_i = angle_init_e + i * an_e
            angle_f = angle_init_e + (i + 1) * an_e
            for r_ in [ri, re]:
                xi = xc - r_ * math.cos(angle_i)
                xf = xc - r_ * math.cos(angle_f)
                zi = r_ * math.sin(angle_i) - zc
                zf = r_ * math.sin(angle_f) - zc
                lines_arch.append({         # Dictionary with the shape of the structure
                    'start': [xi, zi],
                    'end':   [xf, zf],
                    'color': '000000',
                    'width': 0.5,
                })

    lines, vertices = lines_and_points_from_form(forms[0], plot_reactions, cracks, radius, max_width, fix_width, hide_negative=hide_negative)
    total_lines = lines_arch + lines
    linecollection = plotter.draw_lines(total_lines)
    pointcollection = plotter.draw_points(vertices)
    segments = []
    centers = []
    for line in total_lines:
        segments.append([line['start'], line['end']])
    for pt in vertices:
        centers.append(pt['pos'])
    plotter.update_linecollection(linecollection, segments)
    plotter.update_pointcollection(pointcollection, centers, radius=radius)
    plotter.update(pause=delay)
    image = os.path.join(tempfolder, pattern.format(str(img_count)))
    images.append(image)
    plotter.save(image)
    img_count += 1

    for form in forms[1:]:
        lines, vertices = lines_and_points_from_form(form, plot_reactions, cracks, radius, max_width, fix_width, hide_negative=hide_negative)
        total_lines = lines_arch + lines
        plotter.clear_collection(linecollection)
        plotter.clear_collection(pointcollection)
        linecollection = plotter.draw_lines(total_lines)
        pointcollection = plotter.draw_points(vertices)
        segments = []
        centers = []
        for line in total_lines:
            segments.append([line['start'], line['end']])
        for pt in vertices:
            centers.append(pt['pos'])
        plotter.update_linecollection(linecollection, segments)
        plotter.update_pointcollection(pointcollection, centers, radius=radius)
        plotter.update(pause=delay)
        image = os.path.join(tempfolder, pattern.format(str(img_count)))
        images.append(image)
        plotter.save(image)
        img_count += 1

    # plotter.draw_lines(lines)
    # plotter.draw_points(vertices)
    # plotter.draw_points(nodes)

    if save:
        plotter.save_gif(save, images, delay=delay*100)
        # plotter.save(save)
    return plotter


def plot_form_semicirculararch_xz(form, radius=0.05, fix_width=False, max_width=10, simple=False, show_q=True, heights=False, show_edgeuv=False, save=None, thk=0.20, plot_reactions=False, joints=False, cracks=False, yrange = None, linestyle='solid'):
    """ Plot of a 2D diagram in the XZ plane

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    max_width : bool
        Maximum width of the plot.
    simple : bool
        Simple red and blue colour plotting.
    show_q : bool
        Show the force densities on the edges.
    heights : str
        Plot the heights of the nodes.
    show_edgeuv : bool
        Show u,v of the edges.
    thck : float
        Thickness of the structure to plot.
    plot_reactions : bool
        Plot the reaction's extension.
    joints : bool
        Plot joints.
    cracks : bool
        Highlight crack points.
    yrange : list
        If 3D structure sliced, it restricts the yrange to plot.
    linestyle : str
        linestyle for the thrust.

    Returns
    ----------
    obj
        Plotter object.

    """

    i_k = form.index_key()
    gkey_key = form.gkey_key()
    q = [form.edge_attribute((u,v), 'q') for u, v in form.edges_where({'_is_edge': True})]
    qmax = max(abs(array(q)))
    lines = []
    xs = []
    reac_lines = []

    if yrange == None:
        edges_considered = list(form.edges())
        vertices_considered = list(form.vertices())
    else:
        edges_considered = []
        vertices_considered = []
        edges_added = []
        for u, v in form.edges():
            if (yrange[0] <= form.vertex_coordinates(u)[1] <= yrange[1]) and (yrange[0] <= form.vertex_coordinates(v)[1] <= yrange[1]):
                edges_considered.append((u,v))
        for key in form.vertices():
            if yrange[0] < form.vertex_coordinates(key)[1] < yrange[1]:
                vertices_considered.append(key)
                edges_added.append(form.vertex_coordinates(key))
                ngb = list(form.vertex_neighborhood(key))
                print(ngb)
                edges_added.append(form.vertex_coordinates(ngb[1]))
        # print(vertices_considered)
        print('Drawing initial total of vertices:', len(vertices_considered))
        print('Drawing initial total of edges:', len(edges_considered))
        edges_added.sort()
        if len(edges_considered) == 0:
            for i in range(len(edges_added)):

                vertices_considered.append(gkey_key[geometric_key(edges_added[i])])
                if i < len(edges_added)-1:
                    u = gkey_key[geometric_key(edges_added[i])]
                    v = gkey_key[geometric_key(edges_added[i+1])]
                    edges_considered.append((u, v))
        print('Drawing total of edges:', len(edges_considered))
        print(edges_considered)
        print(vertices_considered)
        print('Drawing total of vertices:', len(vertices_considered))


    for key in vertices_considered:
        xs.append(form.vertex_coordinates(key)[0])
        if form.vertex_attribute(key, 'is_fixed') == True:
            x, _, z = form.vertex_coordinates(key)
            if z > 0.0:
                rz = abs(form.vertex_attribute(key, '_rz'))
                rx = form.vertex_attribute(key, '_rx')
                reac_line = [x, z, x + z * rx / rz, 0.0]
                reac_lines.append(reac_line)

    for u, v in edges_considered:
        qi = form.edge_attribute((u, v), 'q')

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
            if form.edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.edge_attribute((u, v), 'is_ind'):
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
            'linestyle': linestyle
        })

    try:
        Re = form.attributes['Re']
        Ri = form.attributes['Ri']
    except:
        Re = 1.1  # (max(xs) - min(xs))/2 + thk/2
        Ri = 0.9  # (max(xs) - min(xs))/2 - thk/2

    xc = (max(xs) - min(xs))/2
    discr = 200
    print('xs:', xc)
    print('Visualisation on Re: {0:.3f} / Ri: {1:.3f}'.format(Re,Ri))

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
                'color': ''.join(['00', '00', '00']),
                'width': max_width,
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
            z = form.vertex_attribute(key, 'ub')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })
        for i in cracks_lb:
            key = i_k[i]
            x, _, _ = form.vertex_coordinates(key)
            z = form.vertex_attribute(key, 'lb')
            vertices.append({
                'pos': [x, z],
                'radius': radius,
                'color': '000000',
            })

    nodes = []
    if radius:
        for key in vertices_considered:
            x, _, z = form.vertex_coordinates(key)
            if form.vertex_attribute(key, 'is_fixed') is True:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '000000',
                    'facecolor': 'aaaaaa',
                })
            if abs(form.vertex_attribute(key, 'ub') - z) < 1e-3:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '008000',
                    'facecolor': '008000',
                })
            if abs(form.vertex_attribute(key, 'lb') - z) < 1e-3:
                nodes.append({
                    'pos': [x, z],
                    'radius': radius,
                    'edgecolor': '0000FF',
                    'facecolor': '0000FF',
                })

    plotter = MeshPlotter(form, figsize=(10, 10))
    # round(form.vertex_attribute(i, 'pz'), 2)
    # if radius:
    #     if heights:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius, text={i: i for i in form.vertices()}) # form.vertex_attribute(i, 'z')
    #     else:
    #         plotter.draw_vertices(facecolor={i: '#aaaaaa' for i in form.vertices_where({'is_fixed': True})},
    #         radius=radius)

    # plotter.draw_vertices(radius= {i : form.vertex_attribute(i, 'px')/100 for i in form.vertices()}) # form.vertex_attribute(i, 'z')

    plotter.draw_lines(lines)
    plotter.draw_points(vertices)
    plotter.draw_points(nodes)

    if save:
        plotter.save(save)

    return plotter


def plot_independents(form, radius=0.05, fix_width=True, width=10, number_ind=True, show_symmetry=False, save=False):
    """ Extended plotting of a FormDiagram focusing on showing independent edges

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    width : bool
        Width of the lines in the plot.
    max_width : float
        Maximum edge width.
    number_ind : bool
        Show or not the numbering on the independent edges.
    show_symmetry : bool
        Show or not the numbering on the symmetrical independent edges.
    save : str
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """

    lines = []
    i = 0

    for u, v in form.edges_where({'_is_edge': True}):
        colour = ['66', '66', '66']
        text = ''
        if form.edge_attribute((u, v), 'is_ind'):
            colour = ['F9', '57', '93']
            if number_ind:
                text = str(i)
                i = i + 1
                if show_symmetry:
                    text = str(form.edge_attribute((u, v), 'sym_key'))

        lines.append({
            'start': form.vertex_coordinates(u),
            'end':   form.vertex_coordinates(v),
            'color': ''.join(colour),
            'width': width,
            'text': text,
        })

    rad_colors = {}
    for key in form.vertices_where({'is_fixed': True}):
        rad_colors[key] = '#aaaaaa'
    for key in form.vertices_where({'rol_x': True}):
        rad_colors[key] = '#ffb733'
    for key in form.vertices_where({'rol_y': True}):
        rad_colors[key] = '#ffb733'

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        plotter.draw_vertices(facecolor=rad_colors, radius=radius)

    # plotter.draw_vertices(keys=[key in form.vertices_where({'is_fixed': True})], radius=10*radius)

    plotter.draw_lines(lines)
    if save:
        plotter.save(save)

    return plotter


def plot_sym_inds(form, radius=0.05, fix_width=True, width=10, save=False):
    """ Extended plotting of a FormDiagram focusing on showing the symmetric relations among independent edges

    Parameters
    ----------
    form : obj
        FormDiagram to plot.
    radius : float
        Radius of vertex markers.
    fix_width : bool
        Fix edge widths as constant.
    width : bool
        Width of the lines in the plot.
    max_width : float
        Maximum edge width.
    save : str
        Path to save the figure, if desired.

    Returns
    ----------
    obj
        Plotter object.

    """

    lines = []
    i = 0

    i_sym_max = 0
    for u, v in form.edges_where({'is_ind': True}):
        i_sym = form.edge_attribute((u, v), 'sym_key')
        if i_sym > i_sym_max:
            i_sym_max = i_sym

    from compas.utilities import rgb_to_hex
    colormap = plt.cm.hsv # gist_ncar nipy_spectral, Set1, Paired coolwarm
    colors = [rgb_to_hex(colormap(i)[:3]) for i in linspace(0, 1.0, i_sym_max + 1)]

    for u, v in form.edges_where({'_is_edge': True}):
        colour = '666666'
        if form.edge_attribute((u, v), 'is_ind'):
            i_sym = form.edge_attribute((u, v), 'sym_key')
            colour = colors[i_sym]

        lines.append({
            'start': form.vertex_coordinates(u),
            'end':   form.vertex_coordinates(v),
            'color': colour,
            'width': width,
            'text': '',
        })

    rad_colors = {}
    for key in form.vertices_where({'is_fixed': True}):
        rad_colors[key] = '#aaaaaa'
    for key in form.vertices_where({'rol_x': True}):
        rad_colors[key] = '#ffb733'
    for key in form.vertices_where({'rol_y': True}):
        rad_colors[key] = '#ffb733'

    plotter = MeshPlotter(form, figsize=(10, 10))
    if radius:
        plotter.draw_vertices(facecolor=rad_colors, radius=radius)

    # plotter.draw_vertices(keys=[key in form.vertices_where({'is_fixed': True})], radius=10*radius)

    plotter.draw_lines(lines)
    if save:
        plotter.save(save)

    return plotter
