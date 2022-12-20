# from compas_tno.diagrams import FormDiagram
# from compas_tno.viewers import Viewer
# from compas_tno.shapes import Shape

# H = 1.0
# L = 2.0
# b = 0.5
# discr = 20

# arch = Shape.create_arch(H=H, L=L, b=b)

# form = FormDiagram.create_arch(H=H, L=L, discretisation=discr)

# view = Viewer(form, arch)
# view.draw_form()
# view.draw_shape()
# view.show()


import compas
from compas.datastructures import Mesh
from compas.geometry import Translation
from compas.geometry import Scale

from compas_view2.app import App

mesh = Mesh.from_obj(compas.get('faces.obj'))
mesh.transform(Translation.from_vector([0.5, 0, 0.1]))

mesh2 = mesh.transformed(Translation.from_vector([-11, 0, 0]))

# =============================================================================
# Visualization
# =============================================================================

viewer = App(width=1600, height=900)
viewer.view.camera.rx = -60
viewer.view.camera.rz = 0
viewer.view.camera.ty = -2
viewer.view.camera.distance = 10

viewer.add(mesh, hide_coplanaredges=False)
viewer.add(mesh2, hide_coplanaredges=True)
viewer.show()
