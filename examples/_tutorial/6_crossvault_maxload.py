from numpy import zeros

from compas.colors import Color
from compas.geometry import distance_point_point_xy
from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.utilities import form_add_lines_support

from compas_viewer import Viewer
import math
from compas.colors import Color
from compas.geometry import Cylinder

# ----------------------------------------
# 0. Vizualization function
# ----------------------------------------

def visualization(vault, form):
    viewer = Viewer()

    # viewer.scene.add(vault.middle, show_lines=False, name="Middle")
    viewer.scene.add(vault.intrados, show_lines=False, name="Intrados", opacity=0.5)
    viewer.scene.add(vault.extrados, show_lines=False, name="Extrados", opacity=0.5)

    edges = list(form.edges_where({"_is_edge": True}))

    max_thick = 0.1
    forces = [form.edge_attribute(edge, "q") * form.edge_length(edge) for edge in edges]
    fmax = math.sqrt(max(abs(max(forces)), abs(min(forces))))

    pipes = []
    for edge in edges:
        qi = form.edge_attribute(edge, "q")
        line = form.edge_line(edge)
        length = line.length
        force = math.sqrt(abs(qi * length))
        radius = force / fmax * max_thick
        pipe = Cylinder.from_line_and_radius(line, radius)
        if force > 1e-3:
            pipes.append(pipe)
        viewer.scene.add(pipe, color=Color.red())

    # viewer.scene.add(pipes, name="Pipes", color=Color.red())

    viewer.show()

# ----------------------------------------
# 1. Shape geometric definition
# ----------------------------------------
spr_angle = 30.0
L = 10.0
thk = 0.50
xy_span = [[0, L], [0, L]]
vault = Shape.create_crossvault(xy_span=xy_span, thk=thk, spr_angle=30)

# ------------------------------------------------------------------------
# 2. Form diagram geometric definition with additional line to supports
# ------------------------------------------------------------------------
discretisation = 8
form = FormDiagram.create_cross_form(xy_span=xy_span, discretisation=discretisation)

load_pos = 3
xc = yc = L / 2
yp = 2.5
yp = yc - load_pos / (discretisation / 2) * yc
for key in form.vertices():
    pt = form.vertex_coordinates(key)
    if distance_point_point_xy(pt, [xc, yp, 0.0]) < 1e-3:
        loaded_node = key
        break

supports = []
for key in form.vertices_where({"is_fixed": True}):
    x, y, z = form.vertex_coordinates(key)
    if y < yc:
        supports.append(key)

print(loaded_node, supports)

form, loaded_node = form_add_lines_support(form, loaded_node=loaded_node, supports=supports)

# plotter = TNOPlotter(form)
# plotter.draw_form(scale_width=False, color=Color.black())
# plotter.draw_supports(color=Color.red())
# plotter.show()

# ------------------------------------------------------------------------
# 4. Define applied load case
# ------------------------------------------------------------------------

n = form.number_of_vertices()
load_direction = zeros((n, 1))
load_direction[loaded_node] = -1.0
print("New Loaded Node:", loaded_node)

# --------------------------------------------
# 5. Maximum load problem and visualisation
# --------------------------------------------
analysis = Analysis.create_max_load_analysis(form, vault, 
                                             load_direction=load_direction, 
                                             max_lambd=300, 
                                             printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

visualization(vault, form)

# viewer = Viewer(form)
# viewer.settings["scale.reactions"] = 0.004
# viewer.draw_thrust()
# viewer.draw_cracks()
# viewer.draw_shape()
# viewer.draw_reactions()

# length = 2.0
# x, y, z = form.vertex_coordinates(loaded_node)
# z += length + 0.1
# arrow = Arrow([x, y, z], [0, 0, -length])
# viewer.app.add(arrow, linecolor=(0, 0, 0), facecolor=(0, 0, 0))

# viewer.show()
