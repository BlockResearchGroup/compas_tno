# ----------------------------------------
# Note: EXAMPLE NEEDS THE INSTALLATION OF THE COMPAS MASONRY VIEWER
# ----------------------------------------

import numpy as np
from compas_masonry.viewers import MasonryViewer

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import distance_point_point_xy
from compas_tna.diagrams import FormDiagram
from compas_tna.envelope import CrossVaultEnvelope
from compas_tno.analysis import Analysis
from compas_tno.utilities import form_add_lines_support

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
form = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)

# ------------------------------------------------------------------------
# 3. Apply Load
# ------------------------------------------------------------------------

load_pos = 4
xc = yc = L / 2
yp = 2.5
yp = yc - load_pos / (n / 2) * yc
for key in form.vertices():
    pt = form.vertex_coordinates(key)
    if distance_point_point_xy(pt, [xc, yp, 0.0]) < 1e-3:
        loaded_node = key
        break

supports = []
for key in form.vertices_where({"is_support": True}):
    coords = form.vertex_coordinates(key)
    if coords is not None:
        x, y, z = coords
        if y < yc:
            supports.append(key)

if load_pos != 0:
    print(loaded_node, supports)
    form, loaded_node = form_add_lines_support(form, loaded_node=loaded_node, supports=supports)

n = form.number_of_vertices()
load_direction = np.zeros((n, 1))
load_direction[loaded_node] = -1.0
print("Loaded Node:", loaded_node)

# --------------------------------------------
# 4. Maximum load problem and visualisation
# --------------------------------------------
analysis = Analysis.create_max_load_analysis(form, vault, load_direction=load_direction, max_lambd=10000, printout=True,)
analysis.apply_selfweight()
analysis.apply_envelope()

analysis.set_up_optimiser()
analysis.run()

viewer = MasonryViewer(formdiagram=form, envelope=vault)
viewer.setup()

# Add load vector to the viewer
x, y, z = form.vertex_coordinates(loaded_node)
load_vector = Vector(0.0, 0.0, -1.0)
anchor_point = Point(x, y, z + 1.0)
viewer.scene.add(load_vector, anchor=anchor_point, name=f"Load Applied_{analysis.optimiser.fopt}", linecolor=Color.black(), linewidth=4.0)

viewer.show()
