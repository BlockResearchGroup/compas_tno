from compas_tno.shapes import Shape

data = {'thk': 0.5, 'discretisation': 10, 't': 0.0, 'xy_span': [[0.0, 10.0], [0.0, 10.0]], 'hc': 5.0, 'hm': None, 'he': None, 'type': 'pointed_crossvault'}

data = {
    'type': 'pointed_crossvault',
    'thk': 0.50,
    'discretisation': 20,
    'xy_span': [[0, 10.0], [0, 10.0]],
    'hc': 5.0,
    'hm': None,
    'he': None,
    'center': [5.0, 5.0],
    'radius': 5.0,
    't': 0.0,
}

shape = Shape.from_library(data)

print(shape)
