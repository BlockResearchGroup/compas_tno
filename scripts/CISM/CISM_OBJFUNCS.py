from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.analysis import Analysis

from compas.geometry import normalize_vector

from numpy import array

import math

discr = 16
xf = 10.0
x0 = 0.0
spr_angle = 30
xc = yc = (x0 + xf)/2
xyspan = [[x0, xf], [x0, xf]]
alpha = 1/math.cos(math.radians(spr_angle))
L = xf * alpha
Ldiff = L - xf
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=discr)
vault = Shape.create_crossvault(xy_span=xyspan_shape, discretisation=discr*2)

vault.ro = 1.0

objective = 't'

if objective == 'min':
    analysis = Analysis.create_minthrust_analysis(form, vault, printout=True)
elif objective == 'max':
    analysis = Analysis.create_maxthrust_analysis(form, vault)
elif objective == 't':
    analysis = Analysis.create_minthk_analysis(form, vault)
elif objective == 'compl_energy':
    vector_supports = []
    for key in form.vertices_where({'is_fixed': True}):
        x, y, z = form.vertex_coordinates(key)
        dXbi = [0, 0, 0]
        if x > xc and y > yc:
            dXbi = normalize_vector([(x - xc), (y - yc), (z - 0.0)])  # outward corner

        vector_supports.append(dXbi)
        form.vertex_attribute(key, 'dXb', dXbi)
    dXb = array(vector_supports)
    analysis = Analysis.create_compl_energy_analysis(form, vault, support_displacement=dXb)

analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
vault0: Shape = Shape.from_formdiagram_and_attributes(form)
analysis.run()

swt = form.lumped_swt()
thrust = form.thrust()
print('T/W:', round(thrust/swt, 4))

swt_0 = vault0.compute_selfweight()
area = vault0.middle.area()
print('Original SWT:', swt_0)
print('Original Area:', area)

view: Viewer = Viewer(form, vault0)
view.draw_thrust(absolute_scale=False)
view.draw_cracks()
view.draw_shape()
view.draw_reactions()
view.show()

# ####### IF WANTS TO REMOVE ELEMENTS

# points = []
# supports = []
# edges = []

# for u, v in form.edges_where({'_is_edge': True}):
#     Xu = form.vertex_coordinates(u)
#     Xv = form.vertex_coordinates(v)

#     if abs(Xu[0] - xc) < 1e-3 and abs(Xv[0] - xc) < 1e-3:
#     # if abs(Xu[1] - Xu[0]) < 1e-3 and abs(Xv[1] - Xv[0]) < 1e-3:
#     # if abs(Xu[1] + Xu[0] - 10.0) < 1e-3 and abs(Xv[1] + Xv[0] - 10) < 1e-3:
#         edges.append((u, v))
#         if u not in points:
#             points.append(u)
#         if v not in points:
#             points.append(v)
#         if form.vertex_attributes(u, 'is_fixed'):
#             supports.append(u)
#         if form.vertex_attributes(v, 'is_fixed'):
#             supports.append(v)

view: Viewer = Viewer(form, shape=vault0)
view.draw_form(edges=edges, cull_negative=True)
# view.draw_cracks(points=points, cull_negative=True)
# view.draw_reactions(supports=supports, extend_reactions=False)
view.draw_shape()
view.show()
