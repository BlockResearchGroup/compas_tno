from compas.geometry import Point
import compas_rhino
import json

# Modify Parameters
titles = ['continuous1', 'continuous2', 'continuous3', 'corners1', 'corners2', 'corners3']
layers = ['Surfaces::' + titles[i] for i in range(len(titles))]
j = 0
for layer in layers:
    data = {'target': {}}
    points = []
    i = 0
    for guid in compas_rhino.get_points(layer=layer):
        point = compas_rhino.rs.PointCoordinates(guid)
        data['target'][i] = list(point)
        i += 1
    print('Found {0} points in layer {1}'.format(i, layer))
    with open('/Users/mricardo/compas_dev/me/freeform/IASS/' + titles[j] + '.json', 'w') as outfile:
        json.dump(data, outfile)
    j += 1