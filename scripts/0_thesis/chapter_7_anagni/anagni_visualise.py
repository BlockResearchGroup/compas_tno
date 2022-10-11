from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from compas_tno.shapes import Shape
from compas.colors import Color

# path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-D3/sangelo_vault_top_final_mesh-D3_with_fill_n.json'
path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-C/sangelo_vault_top_final_mesh-C_with_fill_max.json'
# path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-B5/sangelo_vault_top_final_mesh-B5_with_fill_min.json'

path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-B/sangelo_vault_top_final_mesh-B_with_fill_min.json'
path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-B/sangelo_vault_top_final_mesh-B_with_fill_max.json'

path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-C/sangelo_vault_top_final_mesh-C_with_fill_min.json'
path = '/Users/mricardo/compas_dev/me/anagni/revision/mesh-C/sangelo_vault_top_final_mesh-C_with_fill_max.json'
form = FormDiagram.from_json(path)

# plot = TNOPlotter(form)
# plot.draw_form_independents()
# plot.draw_supports()
# plot.show()

inds = len(list(form.edges_where({'is_ind': True})))
print('NUMBED INDS:', inds)
swt = form.lumped_swt()
thrust = form.thrust()
print(swt, thrust, thrust/swt)

vault = Shape.from_formdiagram_and_attributes(form)

view = Viewer(form, vault)
view.settings['camera.target'] = [2.8, 1.6, 12]
view.settings['camera.distance'] = 18
view.settings['scale.reactions'] = 0.005 * 4
view.settings['scale.loads'] = 0.005 * 4 * 10
view.settings['opacity.shapes'] = 0.3
view.scale_edge_thickness(3.0)
view.draw_form()
view.draw_shape()
view.draw_cracks()
view.show()


view = Viewer(form, vault)
view.settings['camera.target'] = [2.8, 1.6, 12]
view.settings['camera.distance'] = 18
view.settings['scale.reactions'] = 0.005 * 4
view.settings['scale.loads'] = 0.005 * 4 * 10
view.settings['opacity.shapes'] = 0.3
view.settings['color.mesh.middle'] = Color.orange()
view.initiate_app(viewmode="shaded")
# view.scale_edge_thickness(3.0)
# view.draw_form()
view.draw_shape()
view.draw_middle_shape()
# view.draw_cracks()
view.show()
