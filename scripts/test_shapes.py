
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes


# data = {
#     'type': 'crossvault',
#     'thk': 0.5,
#     'density': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

# data = {
#     'type': 'pavillionvault',
#     'thk': 0.5,
#     'density': [50, 50],
#     'xy_span': [[0.0,10.0],[0.0,10.0]],
#     't' : 0.0
# }

data = {
    'type': 'dome',
    'thk': 0.5,
    'density': [50, 50],
    'center': [[5.0,5.0]],
    'radius': 5.0,
    't' : 0.0
}

print( '\nNow B\n')

vault = Shape.from_library(data)

print('\Evaluate Height of some points:')

print('\nMiddle')
print(vault.get_target(5,5))
print(vault.get_target(7.5,7.5))
print(vault.get_target(10.0,5.0))
print(vault.get_target(0.0,5.0))

print('Extrados')
print(vault.get_ub(5,5))
print(vault.get_ub(7.5,7.5))
print(vault.get_ub(10.0,5.0))
print(vault.get_ub(0.0,5.0))

print('Intrados')
print(vault.get_lb(5,5))
print(vault.get_lb(7.5,7.5))
print(vault.get_lb(10.0 - 0.25, 5.0))
print(vault.get_lb(0.0 + 0.25, 5.0))
print(vault.get_lb(10.0,5.0))
print(vault.get_lb(0.0,5.0))

view_shapes(vault).show()
