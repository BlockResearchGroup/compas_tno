from compas_tno.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape

from compas_viewer import Viewer
import math
from compas.colors import Color
from compas.geometry import Cylinder
from compas_tno.viewer import TNOViewer

form = FormDiagram.create_cross_form(discretisation=20)
shape = Shape.create_crossvault(spr_angle=30.0)

analysis = Analysis.create_lp_analysis(form, shape, solver="CVXPY", printout=True)
analysis.apply_selfweight()
analysis.set_up_optimiser()
analysis.run()

TNOViewer(form).show()
