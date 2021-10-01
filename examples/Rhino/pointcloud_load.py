from compas.geometry import Point
import compas_rhino
import json

# Modify Parameters
title = 'sangelo'
layers = ['FormDiagram-Extrados', 'FormDiagram-Intrados']

# Create dictionary and save as .json
points_UB = []
points_LB = []
UBLB = ['UB', 'LB']
data = {UBLB[0]:{}, UBLB[1]:{}}

for j in range(len(layers)):
    i = 0
    for guid in compas_rhino.get_points(layer=layers[j]):
        point = compas_rhino.rs.PointCoordinates(guid)
        data[UBLB[j]][i] = list(point)
        i += 1
    print('Found {0} points in layer {1}'.format(i, layers[j]))


with open('/Users/mricardo/compas_dev/me/anagni/' + title + '.json', 'w') as outfile:
    json.dump(data, outfile)
