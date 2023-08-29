from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter

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
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.optimiser.set_starting_point('current')
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# ----------------------------------------
# 4. Visualise solution
# ----------------------------------------

view = Viewer(form, arch)
view.settings['camera.target'] = [0.5, 0, 0]
view.settings['camera.distance'] = 7.0
view.settings['camera.rx'] = 60
view.draw_form()
view.draw_shape()
view.draw_cracks()
view.show()

# ----------------------------------------
# 5. Create analysis for maximum thrust and visualise
# ----------------------------------------

analysis = Analysis.create_maxthrust_analysis(form, arch)
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.optimiser.set_starting_point('current')
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, arch)
view.settings['camera.target'] = [0.5, 0, 0]
view.settings['camera.distance'] = 7.0
view.settings['camera.rx'] = 60
view.draw_form()
view.draw_shape()
view.draw_cracks()
view.draw_reactions(extend_reactions=True)
view.show()
