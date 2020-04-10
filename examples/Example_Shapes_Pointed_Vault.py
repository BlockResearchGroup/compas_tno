
from compas_tno.shapes import Shape
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_middle


# --------------------------------------------------------------
# -----------PARAMERIC SHAPES (uncomment data) -----------------
# --------------------------------------------------------------

data = {
    'type': 'pointed_crossvault',
    'thk': 0.5,
    'discretisation': [50, 50],
    'xy_span': [[0.0, 6.0], [0.0, 10.0]],
    't': 0.0,
    'hc': 6.0,
    'hm': None,
    'he': None,
}

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
