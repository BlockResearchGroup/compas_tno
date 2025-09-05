import numpy as np
from compas_masonry.viewers import MasonryViewer

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Vector
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
form = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)

# ----------------------------------------
# 3. Define displacement on supports - corner diagonal
# ----------------------------------------

corner_vertices = form.supports()
# Create empty array to add row by row
displacement_array = []
displ_point = [0, 0, 0]
# Go over the fixed vertices and add rows
for vertex in corner_vertices:
    x, y, z = form.vertex_coordinates(vertex)
    if (x - displ_point[0]) ** 2 + (y - displ_point[1]) ** 2 < 1e-3:
        displacement_array.append([-1, -1, 0])
    else:
        displacement_array.append([0, 0, 0])

# Transform to numpy array
displacement_array = np.array(displacement_array)

analysis = Analysis.create_compl_energy_analysis(form, vault, support_displacement=displacement_array, printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

form = analysis.formdiagram
envelope = analysis.envelope

# =============================================================================
# Viz
# =============================================================================

viewer = MasonryViewer(formdiagram=form, envelope=vault)
viewer.setup()

# Add displacement vector to the viewer
for index, vertex in enumerate(form.supports()):
    x, y, z = form.vertex_coordinates(vertex)
    displacement_vector = Vector(displacement_array[index][0], displacement_array[index][1], displacement_array[index][2])
    anchor_point = Point(x, y, z)

    # Add the vector to the viewer with black color and thick line
    if displacement_vector.length > 0.01:
        viewer.scene.add(displacement_vector, anchor=anchor_point, name=f"Displacement_Corner_{vertex}", linecolor=Color.black(), linewidth=4.0)

viewer.show()
