from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

# file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=6.71/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_25.900192083843347.json'

hc_list = [5.00, 5.48, 5.92, 6.32, 6.71, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]
degs = [0, 10, 20, 30, 40]
type_structure = 'pointed_crossvault'
type_formdiagram = 'cross_fd'
discretisation = 20
thk = 0.50
dist = 15.0

for i in range(len(hc_list)):
    for j in range(len(degs)):
        hc = hc_list[i]
        deg = degs[j]
        folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram,'h='+str(hc), 'min_thk')
        if deg:
            folder = os.path.join(folder, 'deg='+str(deg))
        title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
        forms_address = os.path.join(folder, title)
        address_shape = forms_address + '_' + 'shape' + '_thk_' + str(100*thk) + '.json'
        # print(address_shape)
        form = FormDiagram.from_json(address_shape)
        artist = FormArtist(form)

        artist.draw_edges(displacement=[dist*i, dist*j, 0])
        artist.draw_intrados(displacement=[dist*i, dist*j, 0])
        artist.draw_middle(displacement=[dist*i, dist*j, 0])
        artist.draw_extrados(displacement=[dist*i, dist*j, 0])
        artist.draw_extrados(displacement=[dist*i, dist*j, 0])

rs.EnableRedraw(True)
