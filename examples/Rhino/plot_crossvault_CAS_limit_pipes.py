from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

# file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=6.71/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_25.900192083843347.json'

cross_fd_0 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=0/crossvault_cross_fd_discr_14_deg=0_min_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=0/crossvault_cross_fd_discr_14_deg=0_max_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=0/crossvault_cross_fd_discr_14_deg=0_min_thk_40.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=0/crossvault_cross_fd_discr_14_deg=0_max_thk_40.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=0/crossvault_cross_fd_discr_14_deg=0_min_thk_33.46711357370478.json']

cross_fd_20 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=20/crossvault_cross_fd_discr_14_deg=20_min_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=20/crossvault_cross_fd_discr_14_deg=20_max_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=20/crossvault_cross_fd_discr_14_deg=20_min_thk_35.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=20/crossvault_cross_fd_discr_14_deg=20_max_thk_35.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=20/crossvault_cross_fd_discr_14_deg=20_min_thk_0.23078255778614098.json']

cross_fd_40 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=40/crossvault_cross_fd_discr_14_deg=40_min_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=40/crossvault_cross_fd_discr_14_deg=40_max_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=40/crossvault_cross_fd_discr_14_deg=40_max_thk_25.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=40/crossvault_cross_fd_discr_14_deg=40_min_thk_25.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/cross_fd/deg=40/crossvault_cross_fd_discr_14_deg=40_min_thk_8.880742912702008.json']

fan_fd_0 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=0/crossvault_fan_fd_discr_14_deg=0_min_thk_50.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=0/crossvault_fan_fd_discr_14_deg=0_max_thk_50.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=0/crossvault_fan_fd_discr_14_deg=0_min_thk_46.864942983990964.json']

fan_fd_20 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=20/crossvault_fan_fd_discr_14_deg=20_min_thk_50.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=20/crossvault_fan_fd_discr_14_deg=20_max_thk_50.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=20/crossvault_fan_fd_discr_14_deg=20_min_thk_45.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=20/crossvault_fan_fd_discr_14_deg=20_max_thk_45.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=20/crossvault_fan_fd_discr_14_deg=20_min_thk_37.44557859983187.json']

fan_fd_40 = ['/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=40/crossvault_fan_fd_discr_14_deg=40_min_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=40/crossvault_fan_fd_discr_14_deg=40_max_thk_49.99999999999999.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=40/crossvault_fan_fd_discr_14_deg=40_min_thk_30.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=40/crossvault_fan_fd_discr_14_deg=40_max_thk_30.0.json',
'/Users/mricardo/compas_dev/me/shape_comparison/crossvault/fan_fd/deg=40/crossvault_fan_fd_discr_14_deg=40_min_thk_22.218018774419715.json']

amiens = ['/Users/mricardo/compas_dev/me/max_n/amiens_internet/crossvault/fan_fd/min_max/crossvault_fan_fd_discr_14_offset-method_min_thk_28.13275411729205.json']

j0 = 0
dist = 15.0

files_list = cross_fd_0
j = 0.0

#files_list = cross_fd_20
#j = 1.0

#files_list = cross_fd_40
#j = 2.0

#files_list = fan_fd_0
#j = 3.0

#files_list = fan_fd_20
#j = 4.0

#files_list = fan_fd_40
#j = 5.0

files_list = amiens
j = -2.0

area = 0.0


#for i in range(len(files_list)):
for i in range(1):
    file = files_list[-1]
    
    form = FormDiagram.from_json(file)
    form2 = form.copy()
    #thk = form.attributes['thk']
    thk = float('.'.join(file.split('_')[-1].split('.')[:-1]))/100
    artist = FormArtist(form)
    
    artist.color_vertex_extrados = (0, 128, 0)
    artist.radius_sphere = 0.15
    
    swt = 0.0
    area = 0.0
    for key in form.vertices():
        z_middle = form.vertex_attribute(key, 'target')
        swt += form.vertex_attribute(key, 'pz')
        form2.vertex_attribute(key, 'z', z_middle)
    for key in form2.vertices():
        area += form2.vertex_area(key)
    
    print('thk / swt / area:', thk, swt, area)
    factor = (area*0.5)/swt
    print('factor:', factor)
    artist.scale_forces = 0.02 * factor
    artist.draw_edges(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_intrados(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_extrados(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_middle(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_cracks(displacement=[dist*i, dist*j + j0, 0], spheres=True)
    artist.draw_forcepipes(displacement=[dist*i, dist*j + j0, 0])
        
    rs.AddTextDot('{0:.2f}'.format(thk), [-1.0 + dist*i, dist*j + j0 - 1.0, 0.0])

rs.EnableRedraw(True)
