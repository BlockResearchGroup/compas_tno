from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

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
radius = 5.0
thk = 0.50
center = [5, 5]
dome = Shape.create_dome(radius=radius, thk=thk, center=center)
dome.ro = 1.0

# ----------------------------------------
# 2. Form diagram geometric definition
# ----------------------------------------
discretisation = [16, 20]
form = FormDiagram.create_circular_radial_form(radius=radius, discretisation=discretisation)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthk_analysis(form, dome, solver="IPOPT", printout=True)
analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

thk_min = analysis.optimiser.fopt
dome_min = Shape.create_dome(radius=radius, thk=thk_min, center=center)

visualization(dome_min, form)