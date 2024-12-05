import math

from compas_viewer import Viewer

from compas.geometry import Cylinder
from compas.colors import Color
from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

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
discretisation = 10
form = FormDiagram.create_cross_form(xy_span=xy_span, discretisation=discretisation)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

# view = Viewer(form)
# view.show_solution()

# --------------------------------------------
# 4. Maximum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_maxthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

# view = Viewer(form)
# view.show_solution()

# =============================================================================
# Viz
# =============================================================================

viewer = Viewer()

viewer.scene.add(vault.middle, show_lines=False, name="Middle")
viewer.scene.add(vault.intrados, show_lines=False, name="Intrados", opacity=0.5)
viewer.scene.add(vault.extrados, show_lines=False, name="Extrados", opacity=0.5)

edges = list(form.edges_where({"_is_edge": True}))

max_thick = 0.1
forces = [form.edge_attribute(edge, "q") * form.edge_length(edge) for edge in edges]
fmax = math.sqrt(max(abs(max(forces)), abs(min(forces))))

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

viewer.show()
