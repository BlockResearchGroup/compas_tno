from compas.geometry import Point
import compas_rhino
import json

# Modify Parameters
k = 1
layers = ['UB2', 'LB2']

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
        

with open('/Users/mricardo/compas_dev/me/min_thk/pointcloud/nurbs' + str(k) + '.json', 'w') as outfile:
    json.dump(data, outfile)