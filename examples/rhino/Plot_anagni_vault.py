from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.rhino import FormArtist
from compas_tno.rhino import ForceArtist
from compas_rhino.artists import Artist
from compas_rhino.artists import MeshArtist
from compas.geometry import Line
from compas.colors import Color

address = '/Users/mricardo/compas_dev/me/anagni/top_vault_less_discr_form_min.json'
address = '/Users/mricardo/compas_dev/compas_tno/data/form-target.json'
address = '/Users/mricardo/compas_dev/compas_tno/data/form-min.json'

ad_form = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-CISM-2.json'
ad_force = '/Users/mricardo/compas_dev/compas_tno/data/CISM/force-CISM-2.json'

ad_form = '/Users/mricardo/compas_dev/compas_tno/data/form-min.json'

ad_form = '/Users/mricardo/compas_dev/compas_tno/data/CISM/forces/form-after-opt.json'

ad_form = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-temp.json'

ad_form = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-Ecomp-linear-16-thk-0.5-corner-diagonal.json'
ad_form = '/Users/mricardo/compas_dev/compas_tno/data/CISM/form-fan_fd-Ecomp-linear-16-thk-0.5-corner-diagonal.json'

ad_form = '/Users/mricardo/compas_dev/me/anagni/top_vault_less_discr_form_min.json'
ad_form = '/Users/mricardo/compas_dev/me/anagni/revision/top_vault_top_vault_less_discr_min.json'
ad_form = '/Users/mricardo/compas_dev/me/anagni/revision/form_min_good_sol.json'

ad_form = '/Users/mricardo/compas_dev/me/anagni/revision/form_min_good_sol.json'
ad_form = '/Users/mricardo/compas_dev/me/anagni/revision/top_vault_mesh-B_min.json'
ad_form = '/Users/mricardo/compas_dev/me/anagni/revision/sangelo_vault_top_final_mesh-B_min.json'

form = FormDiagram.from_json(ad_form)
force = ForceDiagram.from_json(ad_force)

print(form)
print(force)

#artist = FormArtist(form)
#artist.pipes_scale = 0.0001 * 5
#artist.draw_reactions(scale=0.03)
#artist.draw_forcepipes()
#artist.redraw()

artist = FormArtist(form)
artist.layer = 'Mesh-working_nofill'
#artist.pipes_scale = 0.003
#artist.draw_reactions(scale=0.002)
#artist.pipes_scale = 0.0001
artist.draw_thrust()
#artist.draw_cracks()
#artist.draw_forcepipes()

#artist.draw_thrust()
artist.draw_from_attributes(attribute='ub')
artist.draw_from_attributes(attribute='lb')
artist.redraw()

#scale_loads = 0.1

#for key in form.vertices():
#    pz = form.vertex_attribute(key, 'pz')
#    x, y, z = form.vertex_coordinates(key)
#    z2 = z - scale_loads*pz
#    line = Line([x, y, z], [x, y, z2])
#    artist = Artist(line)
#    artist.draw()


#for u, v in force.edges():
#    print(u, v)
#    Xu = force.vertex_coordinates(u)
#    Xv = force.vertex_coordinates(v)
#    line = Line(Xu, Xv)
#    print(line)
#    artist = Artist(line)
#    artist.draw()

#artist.redraw()

#artist = MeshArtist(force)
#artist.default_color = Color.black()
#artist.draw()
#artist.redraw()


artist.redraw()