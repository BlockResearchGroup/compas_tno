from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

from compas_tno.viewer import TNOViewer

from compas_viewer import Viewer

from compas.colors import Color
from compas.geometry import Cylinder
import math

# ----------------------------------------
# 1. Geometric definition
# ----------------------------------------
H = 1.0
L = 2.0
b = 0.5
discretisation = 20

arch = Shape.create_arch(H=H, L=L, b=b)

# ----------------------------------------
# 2. Form Diagram
# ----------------------------------------

form = FormDiagram.create_arch(H=H, L=L, discretisation=discretisation)

# ----------------------------------------
# 3. Create analysis for minimum thrust
# ----------------------------------------
analysis = Analysis.create_minthrust_analysis(form, arch)
analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
analysis.optimiser.set_starting_point("current")
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# ----------------------------------------
# 4. Visualise solution
# ----------------------------------------

viewer = TNOViewer(form, arch)
viewer.show()

# ----------------------------------------
# 5. Create analysis for maximum thrust and visualise
# ----------------------------------------

analysis = Analysis.create_maxthrust_analysis(form, arch, printout=True)
analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
analysis.optimiser.set_starting_point("current")
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

viewer = TNOViewer(form, arch)
viewer.show()
