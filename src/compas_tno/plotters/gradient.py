from compas_plotters import MeshPlotter
from numpy import array


__all__ = [
    'plot_grad',
]


def plot_grad(form, radius=0.1, fix_width=False, max_width=10, simple=True):
    """ Extended load-path plotting of a the Gradient

    Parameters
    ----------
    form : FormDiagram
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
        qi = form.edge_attribute((u, v), 'dq')

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
            if form.edge_attribute((u, v), 'is_symmetry'):
                colour[1] = 'cc'
            if form.edge_attribute((u, v), 'is_ind'):
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
