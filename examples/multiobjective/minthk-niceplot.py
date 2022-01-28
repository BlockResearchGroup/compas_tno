from compas_tno.diagrams import FormDiagram
import compas_tno
from compas_tno.algorithms import reciprocal_from_form
from compas_tno.shapes.shape import Shape
from compas_tno.viewers import animation_from_optimisation
from compas_tno.viewers import animation_from_section
from compas_tno.viewers import Viewer

# # ---- MINTHK ----

# DATA_FORM = compas_tno.get('form-minthk.json')
# DATA_FORCE = compas_tno.get('force-minthk.json')

# form = FormDiagram.from_json(DATA_FORM)
# shape = Shape.from_library(
#     {
#         'type': 'dome_polar',
#         # 'thk': 0.50,
#         'thk': 0.205,
#         'discretisation': [20, 40],
#         't': 0.0,
#         'center': [5.0, 5.0],
#         'radius': 5.0,
#     })

# pzt = 0
# for key in form.vertices():
#     pz = form.vertex_attribute(key, 'pz')
#     # form.vertex_attribute(key, 'pz', 2 * pz)
#     pzt += form.vertex_attribute(key, 'pz')

# print('Total weight:', pzt)

# qt = 0
# for edge in form.edges():
#     q = form.edge_attribute(edge, 'q')
#     # form.edge_attribute(edge, 'q', q/2)
#     qt += form.edge_attribute(edge, 'q')

# print('Total qs:', qt)

# force = reciprocal_from_form(form)
# force.to_json(DATA_FORCE)
# print('Saved Force @:', DATA_FORCE)

# viewer = Viewer(form, shape=shape)
# SETTINGS = viewer.settings.copy()

# SETTINGS['camera.show.grid'] = False
# SETTINGS['camera.show.axis'] = False
# SETTINGS['camera.distance'] = 20
# SETTINGS['camera.target'] = [0, 0, 5]
# SETTINGS['force.anchor'] = [10, -6, 0]
# SETTINGS['force.scale'] = 0.1
# SETTINGS['force.scale'] = 0.1
# # SETTINGS['scale.reactions'] = 0.05

# viewer.initiate_app()
# viewer.view_thrust()
# viewer.view_cracks()
# viewer.view_force()
# viewer.view_shape()
# viewer.view_reactions()
# viewer.show()


# # # ---- MINTHRUST ----

# DATA_FORM = compas_tno.get('form-minthrust.json')
# DATA_FORCE = compas_tno.get('force-minthrust.json')

# form = FormDiagram.from_json(DATA_FORM)
# shape = Shape.from_library(
#     {
#         'type': 'dome_polar',
#         'thk': 0.50,
#         'discretisation': [20, 40],
#         't': 0.0,
#         'center': [5.0, 5.0],
#         'radius': 5.0,
#     })

# force = reciprocal_from_form(form)
# force.to_json(DATA_FORCE)

# pzt = 0
# for key in form.vertices():
#     pz = form.vertex_attribute(key, 'pz')
#     pzt += form.vertex_attribute(key, 'pz')

# print('Total weight:', pzt)

# qt = 0
# for edge in form.edges():
#     q = form.edge_attribute(edge, 'q')
#     qt += form.edge_attribute(edge, 'q')

# print('Total qs:', qt)

# viewer = Viewer(form, shape=shape)
# SETTINGS = viewer.settings.copy()

# SETTINGS['camera.show.grid'] = False
# SETTINGS['camera.show.axis'] = False
# SETTINGS['camera.distance'] = 20
# SETTINGS['camera.target'] = [0, 0, 5]
# SETTINGS['force.anchor'] = [10, -6, 0]
# SETTINGS['force.scale'] = 0.1
# SETTINGS['force.scale'] = 0.1
# # SETTINGS['scale.reactions'] = 0.05

# viewer.initiate_app()
# viewer.view_thrust()
# viewer.view_cracks()
# viewer.view_force()
# viewer.view_shape()
# viewer.view_reactions()
# viewer.show()


# # # ---- MAXTHRUST ----

# # DATA_FORM = compas_tno.get('form-maxthrust.json')
# # DATA_FORCE = compas_tno.get('force-maxthrust.json')

# # form = FormDiagram.from_json(DATA_FORM)
# # shape = Shape.from_library(
# #     {
# #         'type': 'dome_polar',
# #         'thk': 0.50,
# #         'discretisation': [20, 40],
# #         't': 0.0,
# #         'center': [5.0, 5.0],
# #         'radius': 5.0,
# #     })

# # force = reciprocal_from_form(form)
# # force.to_json(DATA_FORCE)

# # viewer = Viewer(form, shape=shape)
# # SETTINGS = viewer.settings.copy()

# # SETTINGS['camera.show.grid'] = False
# # SETTINGS['camera.show.axis'] = False
# # SETTINGS['camera.distance'] = 20
# # SETTINGS['camera.target'] = [0, 0, 5]
# # SETTINGS['force.anchor'] = [10, -6, 0]
# # SETTINGS['force.scale'] = 0.1
# # SETTINGS['force.scale'] = 0.1
# # # SETTINGS['scale.reactions'] = 0.05

# # viewer.initiate_app()
# # viewer.view_thrust()
# # viewer.view_cracks()
# # viewer.view_force()
# # viewer.view_shape()
# # viewer.view_reactions()
# # viewer.show()


# # ---- COMPL. ENERGY ----

# DATA_FORM = compas_tno.get('form-compl.json')
# DATA_FORCE = compas_tno.get('force-compl.json')

# form = FormDiagram.from_json(DATA_FORM)
# shape = Shape.from_library(
#     {
#         'type': 'dome_polar',
#         'thk': 0.50,
#         # 'thk': 0.205,
#         'discretisation': [20, 40],
#         't': 0.0,
#         'center': [5.0, 5.0],
#         'radius': 5.0,
#     })

# pzt = 0
# for key in form.vertices():
#     pz = form.vertex_attribute(key, 'pz')
#     # form.vertex_attribute(key, 'pz', 2 * pz)
#     pzt += form.vertex_attribute(key, 'pz')

# print('Total weight:', pzt)

# qt = 0
# for edge in form.edges():
#     q = form.edge_attribute(edge, 'q')
#     # form.edge_attribute(edge, 'q', q/2)
#     qt += form.edge_attribute(edge, 'q')

# print('Total qs:', qt)

# force = reciprocal_from_form(form)
# force.to_json(DATA_FORCE)
# print('Saved Force @:', DATA_FORCE)

# viewer = Viewer(form, shape=shape)
# SETTINGS = viewer.settings.copy()

# SETTINGS['camera.show.grid'] = False
# SETTINGS['camera.show.axis'] = False
# SETTINGS['camera.distance'] = 20
# SETTINGS['camera.target'] = [0, 0, 5]
# SETTINGS['force.anchor'] = [10, -6, 0]
# SETTINGS['force.scale'] = 0.1
# SETTINGS['force.scale'] = 0.1
# # SETTINGS['scale.reactions'] = 0.05

# viewer.initiate_app()
# viewer.view_thrust()
# viewer.view_cracks()
# viewer.view_force()
# viewer.view_shape()
# viewer.view_reactions()
# viewer.show()




# ---- GOTHIC MIN THK ----

DATA_FORM = '/Users/mricardo/Documents/ETH/Journals/Benvenuto_2022/key_files/gothic-min_thk-form.json'
DATA_FORCE = '/Users/mricardo/Documents/ETH/Journals/Benvenuto_2022/key_files/gothic-min_thk-force.json'

# from compas_tno.utilities import update_json
# DATA_FORM = update_json(DATA_FORM)
# print('Updated JSON:', DATA_FORM)

form = FormDiagram.from_json(DATA_FORM)
# shape = Shape.from_library(
#     {
#         'type': 'dome_polar',
#         'thk': 0.50,
#         # 'thk': 0.205,
#         'discretisation': [20, 40],
#         't': 0.0,
#         'center': [5.0, 5.0],
#         'radius': 5.0,
#     })

pzt = 0
for key in form.vertices():
    pz = form.vertex_attribute(key, 'pz')
    form.vertex_attribute(key, 'pz', -1 * pz)
    pzt += form.vertex_attribute(key, 'pz')

print('Total weight:', pzt)

qt = 0
for edge in form.edges():
    q = form.edge_attribute(edge, 'q')
    # print(q)
    # form.edge_attribute(edge, 'q', -1 * q)
    qt += form.edge_attribute(edge, 'q')

print('Total qs:', qt)

form.to_json(DATA_FORM)

force = reciprocal_from_form(form)
force.to_json(DATA_FORCE)
print('Saved Force @:', DATA_FORCE)

viewer = Viewer(form, shape=None)
# viewer = Viewer(form, shape=shape)
SETTINGS = viewer.settings

SETTINGS['size.edge.max_thickness'] = 15

SETTINGS['camera.show.grid'] = False
SETTINGS['camera.show.axis'] = False
SETTINGS['camera.distance'] = 20
SETTINGS['camera.target'] = [0, 0, 5]
SETTINGS['force.anchor'] = [10, -6, 0]
SETTINGS['force.scale'] = 0.05
# SETTINGS['scale.reactions'] = 0.05

viewer.initiate_app()
viewer.view_thrust()
viewer.view_cracks()
viewer.view_force()
viewer.view_shape()
# viewer.view_reactions()
viewer.show()
