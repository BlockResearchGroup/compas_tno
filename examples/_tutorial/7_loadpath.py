from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

form = FormDiagram.create_cross_form()
shape = Shape.create_crossvault(spr_angle=30.0)

analysis = Analysis.create_lp_analysis(form, shape, solver='CVXPY', printout=True)
analysis.apply_selfweight()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form)
view.draw_thrust()
view.show()
