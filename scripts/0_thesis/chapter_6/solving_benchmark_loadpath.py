from compas_tno.diagrams import FormDiagram
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.problems import initialize_loadpath
from compas_tno.utilities import apply_selfweight_from_thrust
from compas.colors import Color

form = FormDiagram.create_cross_form(discretisation=20, fix='all')

m = form.number_of_real_edges()
n = form.number_of_supports()

print('Edges, supports', m, n)

apply_selfweight_from_thrust(form, density=20.0)
initialize_loadpath(form, printout=True, solver_convex='MATLAB')
