from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist

path = '/Users/mricardo/compas_dev/me/min_thk/dome/PAPER_CAS/dome_minthk_review.json'
form = FormDiagram.from_json(path)

artist = FormArtist(form)
#artist.draw_thrust()
artist.draw_reactions(scale=0.1)
artist.redraw()