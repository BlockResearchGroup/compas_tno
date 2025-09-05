# ----------------------------------------
# Note: EXAMPLE NEEDS THE INSTALLATION OF THE COMPAS MASONRY VIEWER
# ----------------------------------------

from compas_masonry.viewers import MasonryViewer

from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import DomeEnvelope
from compas_tno.analysis import Analysis

# ----------------------------------------
# 1. Shape geometric definition
# ----------------------------------------

center = (5.0, 5.0)
radius = 5.0
r_oculus = 0.5
thickness = 0.20
min_lb = 0.5
n_hoops = 12
n_parallels = 24

envelope = DomeEnvelope(
    center=center,
    radius=radius,
    thickness=thickness,
    min_lb=min_lb,
    n_hoops=n_hoops * 2,
    n_parallels=n_parallels * 2,
    r_oculus=r_oculus,
)

# ----------------------------------------
# 2. Form diagram geometric definition
# ----------------------------------------

form = FormDiagram.create_circular_radial_spaced(
    center=center,
    radius=radius,
    n_hoops=n_hoops,
    n_parallels=n_parallels,
    r_oculus=r_oculus,
)

for edge in form.edges_on_boundary():
    form.edge_attribute(edge, "_is_edge", False)

# --------------------------------------------
# 3. Minimum thurst solution and visualisation
# --------------------------------------------
analysis = Analysis.create_minthrust_analysis(form, envelope, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

viewer = MasonryViewer(
    formdiagram=form,
    envelope=envelope,
)
viewer.form_max_thk = 0.03
viewer.setup()
viewer.show()
