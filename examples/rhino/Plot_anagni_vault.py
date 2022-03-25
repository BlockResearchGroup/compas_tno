from compas_tno.diagrams import FormDiagram
from compas_tno.rhino import FormArtist

address = '/Users/mricardo/compas_dev/me/anagni/top_vault_less_discr_form_min.json'
form = FormDiagram.from_json(address)

print(form)

artist = FormArtist(form)
artist.pipes_scale = 0.001
artist.draw_reactions(scale=0.03)
artist.draw_forcepipes()
artist.redraw()