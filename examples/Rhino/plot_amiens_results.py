from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

# file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=6.71/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_25.900192083843347.json'
dist = 15.0
j = 0.0

files_list = ['/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_28.13275411729205.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_30.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_32.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_34.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_36.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_38.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_40.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_42.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_max_thk_44.0.json']

files_list2 = ['/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_28.13275411729205.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_30.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_32.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_34.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_36.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_38.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_40.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_42.0.json',
'/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_44.0.json']

j0 = 15.0
#j0= 0.0

files_list = files_list[::-1]
#files_list2 = files_list2[::-1]

for i in range(len(files_list2)):
    file = files_list[i]
    #file = files_list2[i]
    form = FormDiagram.from_json(file)
    artist = FormArtist(form)

    #artist.draw_edges(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_intrados(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_extrados(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_cracks(displacement=[dist*i, dist*j + j0, 0], tol_crack=10e-3)
    #artist.draw_edges(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_forcepipes(displacement=[dist*i, dist*j, 0])
        
    thk = float('.'.join(file.split('_')[-1].split('.')[:-1]))/100
        
    rs.AddTextDot('{0:.2f}'.format(thk), [-1.0 + dist*i, dist*j + j0 - 1.0, 0.0])
        
        

rs.EnableRedraw(True)
