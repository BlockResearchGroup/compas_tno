from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

H = 1.0
L = 2.0
b = 0.5
discretisation = 20

arch = Shape.create_arch(H=H, L=L, b=b)

form = FormDiagram.create_arch(H=H, L=L, discretisation=discretisation)

analysis = Analysis.create_minthrust_analysis(form,
                                              arch,
                                              printout=True)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, arch)
view.draw_form()
view.draw_shape()
view.draw_cracks()
view.show()
