from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.analysis import Analysis


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
analysis = Analysis.create_minthk_analysis(form, dome, solver='IPOPT', printout=True)
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

thk_min = analysis.optimiser.fopt
dome_min = Shape.create_dome(radius=radius, thk=thk_min, center=center)

view = Viewer(form, dome_min)
view.scale_edge_thickness(5.0)
view.draw_form()
view.draw_shape()
view.draw_reactions(extend_reactions=True)
view.draw_cracks()
view.show()
