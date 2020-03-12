
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.viewers.shapes import view_intrados

# --------------------------------------------------------------
# -----------PARAMERIC SHAPES (uncomment data) -----------------
# --------------------------------------------------------------

# WIP ADD CIRCULAR ARCH

# data = {
#     'type': 'crossvault',
#     'thk': 0.5,
#     'discretisation': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

# data = {
#     'type': 'pavillionvault',
#     'thk': 0.5,
#     'discretisation': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

# data = {
#     'type': 'dome',
#     'thk': 0.5,
#     'discretisation': [50, 50],
#     'center': [5.0, 5.0],
#     'radius': 5.0,
#     't' : 0.0
# }

data = {
    'type': 'arch',
    'H': 0.50,
    'L': 2.0,
    'thk': 0.10,
    'discretisation': 13,
    'b': 0.5,
    't': 0.4,
}

vault = Shape.from_library(data)

print('\Evaluate Height of some points:')

points = [
    [5, 5],
    [7.5, 7.5],
    [10.0, 5.0],
    [0.0, 5.0],
]

points = [
    [0,0],
    [1.0,0.0],
    [2.0,0.0],
    [0.5,0.0],
]

for pt in points:
    print('Point:', pt, 'evaluated on Target / Extrados / Intrados:', pt)
    print(vault.get_middle(pt[0], pt[1]))
    print(vault.get_ub(pt[0], pt[1]))
    print(vault.get_lb(pt[0], pt[1]))

view_shapes(vault).show()
view_intrados(vault).show()
