from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.analysis import Analysis
from compas.geometry import Vector, Point
from numpy import array


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

# ----------------------------------------
# 3. Define displacement field
# ----------------------------------------

vector_supports = []
vectors_plot = []
base_plot = []
xc = center[0]

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    dXbi = [0, 0, 0]
    if x - xc > 0.1:
        dXbi = [1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))
    if x - xc < -0.1:
        dXbi = [-1, 0, 0]
        vectors_plot.append(Vector(*dXbi))
        base_plot.append(Point(x, y, z))

    vector_supports.append(dXbi)

dXb = array(vector_supports)

# --------------------------------------------
# 4. Create analysis, run and visualise
# --------------------------------------------
analysis = Analysis.create_compl_energy_analysis(form, dome, solver='IPOPT', support_displacement=dXb, printout=True)
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, dome)
view.scale_edge_thickness(5.0)
view.draw_form()
view.draw_shape()
view.draw_reactions(extend_reactions=True)
view.draw_cracks()
for i in range(len(vectors_plot)):
    vector = vectors_plot[i]
    base = base_plot[i]
    view.draw_vector(vector=vector, base=base)
view.show()
