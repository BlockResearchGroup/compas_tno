
from compas_tno.shapes import Shape
from compas_tno.viewers.shapes import view_shapes
import math

# --------------------------------------------------------------
# -----------PARAMERIC SHAPES (uncomment data) -----------------
# --------------------------------------------------------------

# data = {
#     'type': 'crossvault',
#     'thk': 0.5,
#     'discretisation': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

# data = {
#     'type': 'pointed_crossvault',
#     'thk': 0.5,
#     'discretisation': [50, 50],
#     'xy_span': [[0.0, 10.0], [0.0, 10.0]],
#     't': 0.0,
#     'hc': 5.0*math.sqrt(2),
#     'hm': None,
#     'he': None,
# }

# data = {
#     'type': 'pavillionvault',
#     'thk': 0.20,
#     'discretisation': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

data = {
    'type': 'pointed_crossvault',
    'thk': 0.5,
    'discretisation': [50, 50],
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0,
    'hc': 5.00,
    'hm': [8.66]*4,
    'he': [5.0, 5.0, 5.0, 5.0],
}

# k1 = 0.0
# k2 = 0.5  # Percentage of the dome to consider
# data = {
#     'type': 'dome_spr',
#     'thk': 0.25,
#     'discretisation': [25, 50],
#     'center': [5.0, 5.0],
#     'radius': 5.0,
#     'theta': [k1 * pi/2, k2 * pi/2],
#     't': 0.0,
# }

# data = {
#     'type': 'dome',
#     'thk': 0.15,
#     'discretisation': [50, 50],
#     't' : 0.0,
#     'center': [5.0, 5.0],
#     'radius': 5.0,
# }

# data = {
#     'type': 'domicalvault',
#     'xy_span': [[0.0, 10.0], [0.0, 10.0]],
#     'thk': 0.50,
#     'discretisation': [50, 50],
#     # 'center': [5.0, 5.0],
#     # 'radius': 8.0,
#     't' : 0.0
# }

# k1 = 0.0
# k2 = 0.5  # Percentage of the dome to consider

# data = {
#     'type': 'dome_spr',
#     'thk': 0.25,
#     'discretisation': [25, 50],
#     'center': [5.0, 5.0],
#     'radius': 5.0,
#     'theta': [k1 * math.pi/2, k2 * math.pi/2],
#     't': 0.0,
# }

# data = {
#     'type': 'arch',
#     'H': 0.50,
#     'L': 2.0,
#     'thk': 0.10,
#     'discretisation': 13,
#     'b': 0.5,
#     't': 0.4,
#     'x0': 0.0,
# }

vault = Shape.from_library(data)

print('\Evaluate Height of some points:')

points = [
    [5, 5],
    [7.5, 7.5],
    [10.0, 5.0],
    [0.0, 5.0],
]

for pt in points:
    print('Point:', pt, 'evaluated on Target / Extrados / Intrados:', pt)
    print(vault.get_middle(pt[0], pt[1]))
    print(vault.get_ub(pt[0], pt[1]))
    print(vault.get_lb(pt[0], pt[1]))

view_shapes(vault).show()
