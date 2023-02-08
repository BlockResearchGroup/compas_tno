from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.analysis import Analysis


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

view = Viewer(form)
view.settings['scale.reactions'] = 0.004
view.show_solution()

# --------------------------------------------
# 4. Maximum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_maxthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form)
view.settings['scale.reactions'] = 0.004
view.show_solution()
