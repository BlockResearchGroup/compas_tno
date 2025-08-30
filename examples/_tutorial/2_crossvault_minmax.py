import math

from compas_viewer import Viewer

from compas.colors import Color
from compas.geometry import Cylinder
from compas.geometry import Line
from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewer import TNOViewer

# ----------------------------------------
# 1. Shape geometric definition
# ----------------------------------------
spr_angle = 30.0
L = 10.0
thk = 0.50
xy_span = [[0, L], [0, L]]
vault = Shape.create_crossvault(xy_span=xy_span, thk=thk, spr_angle=30)

# ----------------------------------------
# 2. Form diagram geometric definition
# ----------------------------------------
discretisation = 14
form = FormDiagram.create_cross_form(xy_span=xy_span, discretisation=discretisation)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

TNOViewer(form, vault).show()

# # --------------------------------------------
# # 4. Maximum thurst solution and visualisation
# # --------------------------------------------
# analysis = Analysis.create_maxthrust_analysis(form, vault, printout=True)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.set_up_optimiser()
# analysis.run()

# TNOViewer(form, vault).show()


viewer = Viewer()

viewer.scene.add(vault.middle, show_lines=False, name="Middle")
viewer.scene.add(vault.intrados, show_lines=False, name="Intrados", opacity=0.5)
viewer.scene.add(vault.extrados, show_lines=False, name="Extrados", opacity=0.5)

edges = list(form.edges_where({"_is_edge": True}))

max_thick = 0.05
forces = [form.edge_attribute(edge, "q") * form.edge_length(edge) for edge in edges]
fmax = math.sqrt(max(abs(max(forces)), abs(min(forces))))

# force pipes

pipes = []
for edge in edges:
    q = form.edge_attribute(edge, "q")
    line = form.edge_line(edge)
    length = line.length
    force = math.sqrt(abs(q * length))
    radius = force / fmax * max_thick
    pipe = Cylinder.from_line_and_radius(line, radius)
    if force > 1e-3:
        pipes.append(pipe)

viewer.scene.add(pipes, name="Pipes", color=Color.red())

# hinge points

top = []
above = []

bottom = []
below = []

for vertex in form.vertices():
    lb = form.vertex_attribute(vertex, "lb")
    ub = form.vertex_attribute(vertex, "ub")
    point = form.vertex_point(vertex)

    if abs(ub - point.z) < 1e-3:
        top.append(point)
    elif abs(lb - point.z) < 1e-3:
        bottom.append(point)

    if point.z > ub:
        above.append(point)
    elif point.z < lb:
        below.append(point)

viewer.scene.add(top, name="Extrados cracks", pointcolor=Color.cyan(), pointsize=20)
viewer.scene.add(above, name="Above", pointcolor=Color.green(), pointsize=20)

viewer.scene.add(bottom, name="Intrados crack", pointcolor=Color.blue(), pointsize=20)
viewer.scene.add(below, name="Below", pointcolor=Color.green(), pointsize=20)

# reactions

lines = []

scale = 1e-2

for vertex in form.vertices_where({"is_fixed": True}):
    point = form.vertex_coordinates(vertex)
    rx = form.vertex_attribute(vertex, "_rx") * scale
    ry = form.vertex_attribute(vertex, "_ry") * scale
    rz = form.vertex_attribute(vertex, "_rz") * scale

    line = Line.from_point_and_vector(point, [-rx, -ry, -rz])
    lines.append(line)

viewer.scene.add(lines, name="Reactions", linecolor=Color.green(), linewidth=5)

viewer.show()
