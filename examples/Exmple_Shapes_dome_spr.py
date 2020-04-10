
from compas_tno.shapes import Shape
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_middle
from math import pi

# --------------------------------------------------------------
# -----------PARAMERIC SHAPES (uncomment data) -----------------
# --------------------------------------------------------------

k1 = 0.0
k2 = 0.5  # Percentage of the dome to consider

data = {
    'type': 'dome_spr',
    'thk': 0.25,
    'discretisation': [25, 50],
    'center': [5.0, 5.0],
    'radius': 5.0,
    'theta': [k1 * pi/2, k2 * pi/2],
    't': 0.0,
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
