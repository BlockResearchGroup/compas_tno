from compas.geometry import Point
import compas_rhino
import json

# Modify Parameters
title = 'jeronimos2'
layers = ['Pts - Copy']

# Create dictionary and save as .json
points_target = []
data = {'target':{}}

for j in range(len(layers)):
    i = 0
    for guid in compas_rhino.get_points(layer=layers[j]):
        point = compas_rhino.rs.PointCoordinates(guid)
        data['target'][i] = list(point)
        i += 1
        

with open('/Users/mricardo/compas_dev/me/min_thk/pointcloud/' + title + '.json', 'w') as outfile:
    json.dump(data, outfile)