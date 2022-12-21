from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.problems import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_shape
from compas.colors import Color
import math

spr_angle = 30.0
R_over_L = 0.5
discr = 16
x0, xf = 0, 10.0
thk = 0.5

form = FormDiagram.create_cross_form(discretisation=discr)

cos = math.cos(math.radians(spr_angle))
fact = 2 * (R_over_L * (cos - 1) + 1/2)
L = (xf - x0)/fact
R = L * R_over_L
Ldiff = L - (xf - x0)
xyspan_shape = [[-Ldiff/2 + x0, xf + Ldiff/2], [-Ldiff/2 + x0, xf + Ldiff/2]]

hc = math.sqrt(R**2 - (R - L/2)**2)

vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk)
vault.ro = 0.1  # really important to IPOPT solution

analysis = Analysis.create_minthrust_analysis(form, vault, printout=True, solver='IPOPT')
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

view = Viewer(form)
view.show_solution()

pzt = form.lumped_swt()
T = form.thrust()

thrust_over_weight = T/pzt

print('T/W:', thrust_over_weight)
