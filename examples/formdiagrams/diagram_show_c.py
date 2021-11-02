from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form


def plot_c_bounds(form, clist=[0.1, 0.25, 0.50]):
    for c in clist:
        lines = []

        for key in form.vertices_where({'is_fixed': False}):
            x, y, _ = form.vertex_coordinates(key)
            print(x, y, c)
            pta = [x - c, y - c, 0.0]
            ptb = [x + c, y - c, 0.0]
            ptc = [x + c, y + c, 0.0]
            ptd = [x - c, y + c, 0.0]
            linesbounds = [[pta, ptb],
                       [ptb, ptc],
                       [ptc, ptd],
                       [ptd, pta],
                      ]
            for line in linesbounds:
                lines.append({
                'start': line[0],
                'end':   line[1],
                'color': '333333',
            })


        from compas_plotters import MeshPlotter
        plotter = MeshPlotter(form)
        plotter.draw_edges()
        plotter.draw_lines(lines)
        plotter.draw_vertices(keys=[key for key in form.fixed()], facecolor=(0, 0, 0))
        plotter.show()

# ------------------------------------------------
# --------- CREATING CROSS FORM DIAGRAM ----------
# ------------------------------------------------

data = {
    'type': 'cross_fd',
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)

plot_c_bounds(form, clist = [0.1, 0.25, 0.50])


# ------------------------------------------------
# --------- CREATING FAN FORM DIAGRAM ------------
# ------------------------------------------------

data = {
    'type': 'fan_fd',
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': [10, 10],
    'fix': 'corners',
}

form = FormDiagram.from_library(data)

plot_c_bounds(form, clist = [0.1, 0.25, 0.50])
