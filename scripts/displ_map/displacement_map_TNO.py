
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities.form import displacement_map_parabola
from compas_tno.utilities.form import displacement_map_4parabolas
from compas_tno.utilities.form import displacement_map_paraboloid
from compas_tno.problems import initialise_form
from compas.colors import Color

from scipy.sparse import diags
from numpy.linalg import matrix_rank
from numpy import vstack
from numpy import hstack
from compas.geometry import Vector
from compas.geometry import Point

form = FormDiagram.create_cross_form()

M = initialise_form(form)

# dX = displacement_map_parabola(form)
dX = displacement_map_4parabolas(form, tol=0.1)
# dX = displacement_map_paraboloid(form, radius=10.0)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False, color=Color.black())

delta = 0.5
for i, key in enumerate(form.vertices()):
    x, y, z = form.vertex_coordinates(key)
    dx, dy = dX[i]
    r = Vector(delta*dx, delta*dy, 0.0)
    pt = Point(x, y, z)
    plotter.app.add(r, point=pt, color=Color.blue())

plotter.show()

Ud = diags(M.C @ dX[:, 0])
Vd = diags(M.C @ dX[:, 1])

Edx = M.Cit @ Ud
Edy = M.Cit @ Vd

Ed = vstack([Edx.todense(), Edy.todense()])

print('Ed shape', Ed.shape)
print('E shape', M.E.shape)

Edr = matrix_rank(Ed)
Er = matrix_rank(M.E)

print('Rank Edr and Er is:', Edr, Er)

E_Edr = hstack([M.E, Ed])
E_Edr2 = vstack([M.E, Ed])

rank_together = matrix_rank(E_Edr)
rank_together2 = matrix_rank(E_Edr2)

print('Rank together is:', rank_together)
print('Rank together is:', rank_together2)

for i, key in enumerate(form.vertices()):
    x, y, z = form.vertex_coordinates(key)
    x = x + dX[i, 0] * delta
    y = y + dX[i, 1] * delta
    form.vertex_attributes(key, 'xy', [x, y])


plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False, color=Color.black())
plotter.show()
