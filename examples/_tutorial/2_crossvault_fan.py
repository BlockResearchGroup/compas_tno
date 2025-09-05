# ----------------------------------------
# Note: EXAMPLE NEEDS THE INSTALLATION OF THE COMPAS MASONRY VIEWER
# ----------------------------------------

from compas_masonry.viewers import MasonryViewer

from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import CrossVaultEnvelope
from compas_tno.analysis import Analysis

# ----------------------------------------
# 1. Shape geometric definition
# ----------------------------------------
L = 10.0
thk = 0.50
x_span = (0, L)
y_span = (0, L)
vault = CrossVaultEnvelope(x_span=x_span, y_span=y_span, thickness=thk)

# ----------------------------------------
# 2. Form diagram geometric definition
# ----------------------------------------
n = 14
form = FormDiagram.create_fan(x_span=x_span, y_span=y_span, n_fans=n, n_hoops=n)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

viewer = MasonryViewer(formdiagram=form, envelope=vault)
viewer.setup()
viewer.show()

# --------------------------------------------
# 4. Maximum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_maxthrust_analysis(form, vault, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

viewer = MasonryViewer(formdiagram=form, envelope=vault)
viewer.setup()
viewer.show()
