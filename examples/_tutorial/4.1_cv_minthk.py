from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

from compas_tno.viewer import TNOViewer

# ----------------------------------------
# 1. Shape geometric definition
# ----------------------------------------
spr_angle = 30.0
L = 10.0
thk = 0.20
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
# analysis = Analysis.create_minthk_analysis(form, vault, solver="IPOPT", printout=True)
analysis = Analysis.create_minthk_analysis(form, vault, solver="SLSQP", printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

print("Minimum thickness:", analysis.optimiser.fopt)

thin_shape = Shape.from_formdiagram_and_attributes(form)
TNOViewer(form, thin_shape).show()