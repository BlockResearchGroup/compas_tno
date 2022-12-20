from compas_rhino.artists import Artist
import compas_rhino
import json


folder_mod = '/Users/mricardo/compas_dev/me/max_load/dome/quadspan/modified_diagrams/'
folder_mod = '/Users/mricardo/compas_dev/me/images/'

json_file = folder_mod + 'solution_D.json'
json_file = '/Users/mricardo/compas_dev/me/minmax/dome/thesis/dome_[20, 16]_thk_0.5_max_raw.json'
json_file = folder_mod + 'fan_lp_raw.json'

f = open(json_file)
data = json.load(f)

lines = data['lines']
points = data['points']

compas_rhino.draw_lines(lines)
compas_rhino.draw_points(points)
     
compas_rhino.redraw()

