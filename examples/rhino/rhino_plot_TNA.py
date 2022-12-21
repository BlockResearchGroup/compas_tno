from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs


path = '/Users/mricardo/compas_dev/me/images/'
file_json = path + 'cross_lp_0.json'
file_raw = path + 'fan_raw.json'

form = FormDiagram.from_json(file_json)
artist = FormArtist(form)

#artist.draw_from_attributes('ub')
#artist.draw_from_attributes('lb')

#artist.draw_mesh()
#artist.draw_reactions(scale=0.3)

rs.EnableRedraw(True)

print(a)

for i in [0, 2, 6]:

    file = '/Users/mricardo/compas_dev/compas_tno/data/form_q={}.json'.format(i)

    form = FormDiagram.from_json(file)
    artist = FormArtist(form)

    artist.pipes_scale = 0.001

    artist.draw_forcepipes(tol=1e-4)
    artist.draw_reactions(draw_as_pipes=False)
    #artist.draw_cracks()
    artist.draw_thrust()

    rs.EnableRedraw(True)