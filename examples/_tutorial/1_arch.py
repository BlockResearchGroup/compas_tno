from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewer import TNOViewer

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
analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
analysis.optimiser.set_starting_point("current")
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

viewer = TNOViewer(form, arch)
viewer.show()


# # ----------------------------------------
# # 4. Create analysis for maximum thrust and visualise
# # ----------------------------------------

# analysis = Analysis.create_maxthrust_analysis(form, arch, printout=True)
# analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
# analysis.optimiser.set_starting_point("current")
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()
# analysis.run()

# viewer = TNOViewer(form, arch)
# viewer.show()

# ----------------------------------------
# 5. Create minimum thickness analysis and visualise
# ----------------------------------------

analysis = Analysis.create_minthk_analysis(form, arch, printout=True)
analysis.optimiser.set_constraints(["funicular", "envelope", "reac_bounds"])
analysis.optimiser.set_starting_point("current")
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

minarch = Shape.create_arch(H=H, L=L, b=b, thk=analysis.optimiser.fopt, t=0.0)

viewer = TNOViewer(form, minarch)
viewer.settings["form_max_thk"] = 0.02
viewer.show()