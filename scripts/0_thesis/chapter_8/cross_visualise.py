import compas_tno
import math
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities import apply_envelope_from_shape

discretisation = 14
thk = 0.5
span = 10.0
spr_angle = 30.0
alpha = 1/math.cos(math.radians(spr_angle))
L = span * alpha
Ldiff = L - span
xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]
shape = Shape.create_crossvault(xy_span=xyspan_shape, thk=thk)

for slide in [0.5]:
# for slide in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    path = '/Users/mricardo/compas_dev/me/hor-loads/cross/slide_diagram/cross_fd_discr_14_tappered_True_slide_{}_lambdh.json'.format(slide)
    path = compas_tno.open_dialog()
    form = FormDiagram.from_json(path)

    W = form.lumped_swt()
    print('SWT:', W)

    for key in form.fixed():
        rx = form.vertex_attribute(key, '_rx')
        ry = form.vertex_attribute(key, '_ry')
        rz = form.vertex_attribute(key, '_rz')
        x, y, z = form.vertex_coordinates(key)
        print(x, y, rx/W, ry/W, rz/W)

    # plot = TNOPlotter(form)
    # plot.draw_form_independents()
    # plot.draw_supports()
    # plot.show()

    # lambdh = form.attributes['lambdh']
    # print('Horizontal loads:', lambdh)
    # print(slide, lambdh)

    base = FormDiagram.create_cross_form(discretisation=50)
    apply_envelope_from_shape(base, shape)
    shape_fixed = Shape.from_formdiagram_and_attributes(base)
    # shape_fixed = Shape.from_formdiagram_and_attributes(form)

    view = Viewer(form, shape_fixed)
    view.scale_edge_thickness(10.0)
    view.draw_thrust()
    view.draw_reactions()
    view.draw_shape()
    # view.draw_reactions(scale=0.5)
    view.draw_cracks()
    view.show()
