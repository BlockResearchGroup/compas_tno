import math

from compas_viewer import Viewer

from compas.colors import Color
from compas.geometry import Cylinder
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
form = FormDiagram.create_fan_form(xy_span=xy_span, discretisation=discretisation)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

TNOViewer(form, vault).show()

# --------------------------------------------
# 4. Maximum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_maxthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

TNOViewer(form, vault).show()