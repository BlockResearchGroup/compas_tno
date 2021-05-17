from compas_tna.diagrams import FormDiagram
from compas_tno.rhino import FormArtist
import rhinoscriptsyntax as rs
import os

# file = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=6.71/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_25.900192083843347.json'

files_fan = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=11.4027/min_thk/pointed_crossvault_fan_fd_discr_14_min_thk_32.59454210856643.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=9.7364/min_thk/deg=10/pointed_crossvault_fan_fd_discr_14_min_thk_28.040144250518782.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=7.9422/min_thk/deg=20/pointed_crossvault_fan_fd_discr_14_min_thk_21.06156962131167.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=6.8896/min_thk/deg=30/pointed_crossvault_fan_fd_discr_14_min_thk_15.38728690454782.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/fan_fd/R=6.2453/min_thk/deg=40/pointed_crossvault_fan_fd_discr_14_min_thk_11.060030746125904.json']
files_cross = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=7.2299/min_thk/pointed_crossvault_cross_fd_discr_14_min_thk_24.188374271079656.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.79/min_thk/deg=10/pointed_crossvault_cross_fd_discr_14_min_thk_20.339922896476374.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.1147/min_thk/deg=20/pointed_crossvault_cross_fd_discr_14_min_thk_13.858171219414173.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.6437/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_8.540529922084234.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.3306/min_thk/deg=40/pointed_crossvault_cross_fd_discr_14_min_thk_4.901892767206777.json']
files_crossbraced = ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/topology-crossbraced/R=6.5/min_thk/deg=20/pointed_crossvault_topology-crossbraced_discr_14_min_thk_15.22778323411051.json']

dome = ['/Users/mricardo/compas_dev/me/min_thk/dome/radial_fd/dome_radial_fd_discr_[20, 16]_min_thk_t_0.20454691527171837.json']
dome = ['/Users/mricardo/compas_dev/me/min_thk/dome/radial_fd/min_max/dome_radial_fd_discr_[20, 16]qmax=5.json']


# ['/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=7.2299/min_thk/pointed_crossvault_cross_fd_discr_14_min_thk_24.188374782556156.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.79/min_thk/deg=10/pointed_crossvault_cross_fd_discr_14_min_thk_20.339922896476374.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=6.1147/min_thk/deg=20/pointed_crossvault_cross_fd_discr_14_min_thk_13.878393453659262.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.6437/min_thk/deg=30/pointed_crossvault_cross_fd_discr_14_min_thk_8.753882727013437.json', '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/R=5.3306/min_thk/deg=40/pointed_crossvault_cross_fd_discr_14_min_thk_4.928791629011682.json']

j0 = -25.0
dist = 25.0

files_list = files_fan
j = 0.0

files_list = files_cross
j = 1.0

files_list = files_crossbraced
j = 2.0

files_list = dome
j = 0.0

#files_list = fan_fd_0
#j = 3.0

#files_list = fan_fd_20
#j = 4.0

#files_list = fan_fd_40
#j = 5.0

area = 0.0


for i in range(len(files_list)):
    file = files_list[i]
    form = FormDiagram.from_json(file)
    form2 = form.copy()
    # thk = form.attributes['thk']
    thk = 0.5  # BY HAND
    # thk = float('.'.join(file.split('_')[-1].split('.')[:-1]))/100
    print('thickness', thk)
    artist = FormArtist(form)
    artist.color_vertex_extrados = (0, 128, 0)
    artist.radius_sphere = 0.15
    swt = 0.0
    area = 0.0
    for key in form.vertices():
        try:
            z_middle = form.vertex_attribute(key, 'target')[0]
        except:
            z_middle = form.vertex_attribute(key, 'target')
        swt += form.vertex_attribute(key, 'pz')
        form2.vertex_attribute(key, 'z', z_middle)
    for key in form2.vertices():
        area += form2.vertex_area(key)
    
    print('thk / swt / area:', thk, swt, area)
    factor = (area*0.5)/swt
    print('factor:', factor)
    
    tol_cracks = 0.01 * thk/100
    
    artist.scale_forces = 0.02 * factor
    artist.draw_edges(displacement=[dist*i, dist*j + j0, 0])
    artist.draw_reactions(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_intrados(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_extrados(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_middle(displacement=[dist*i, dist*j + j0, 0])
    #artist.draw_cracks(displacement=[dist*i, dist*j + j0, 0], spheres=True, tol_cracks=tol_cracks)
    #artist.draw_forcepipes(displacement=[dist*i, dist*j + j0, 0])
    
    # make a draw reaction function
    
    rs.AddTextDot('{0:.2f}'.format(thk), [-1.0 + dist*i, dist*j + j0 - 1.0, 0.0])

rs.EnableRedraw(True)
