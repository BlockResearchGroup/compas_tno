from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_bounds_reactions

paths = [
    # '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_compression_neg.json',
    '/Users/mricardo/compas_dev/me/minmax/dome/thesis/dome_[20, 16]_thk_0.5_min.json',
    # '/Users/mricardo/compas_dev/me/minmax/dome/thesis/dome_[20, 16]_thk_0.5_max.json'
]

for path in paths:

    form = FormDiagram.from_json(path)
    thk = form.attributes['thk']

    shape = Shape.create_dome(thk=thk)

    W = form.lumped_swt()
    print(' Weight:', W)

    apply_envelope_from_shape(form, shape)
    apply_bounds_reactions(form, shape)

    view = Viewer(form, shape=shape)
    view.scale_edge_thickness(7.29)
    view.settings['camera.show.grid'] = False
    view.settings['camera.distance'] = 35
    view.settings['camera.target'] = [5, 5, 0]
    view.settings['camera.rz'] = 45
    view.settings['camera.rx'] = 60
    view.draw_form(absolute_scale=True, cull_negative=True)
    view.draw_cracks(cull_negative=True)
    view.draw_reactions(extend_reactions=True)
    view.draw_shape()
    view.show()

    data = view.to_objects()

    print(data)


