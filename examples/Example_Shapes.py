
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes

# --------------------------------------------------------------
# -----------PARAMERIC SHAPES (uncomment data) -----------------
# --------------------------------------------------------------


data = {
    'type': 'crossvault',
    'thk': 0.5,
    'discretisation': [50, 50],
    'xy_span': [[0.0,8.0],[0.0,10.0]],
    't' : 0.0
}

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

vault = Shape.from_library(data)

print('\Evaluate Height of some points:')

points = [
    [5,5],
    [7.5,7.5],
    [10.0,5.0],
    [0.0,5.0],
]

for pt in points:
    print('Point:', pt, 'evaluated on Target / Extrados / Intrados:', pt)
    print(vault.get_middle(pt[0],pt[1]))
    print(vault.get_ub(pt[0],pt[1]))
    print(vault.get_lb(pt[0],pt[1]))

view_shapes(vault).show()
