from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
# from compas_tno.diagrams.form import overview_forces
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json


# Plot Thrust Network
i = 0
# shapes = [2,4,5,7,8,10,11]
# shapes = range(2,9)

for j in [2]: # j = 1

    if j == 0 or j == 1:
        shapes = [2,4,5,7,8,10,11]
    else:
        shapes = [8] #range(2,8)

    i = 1

    for pat in ['c']:

        # fnm = '/Users/mricardo/compas_dev/me/bestfit/crossvault/discretize/0'+str(j)+'_0'+str(i)+'_fit_crossvault.json'
        # fnm = '/Users/mricardo/compas_dev/me/bestfit/sixpartite/discretize/0'+str(j)+'_0'+str(i)+'_fit_sixpartite.json'

        # Load

        # file = '/Users/mricardo/compas_dev/me/minmax/2D_Arch/01.json'
        # fnm = '/Users/mricardo/compas_dev/me/loadpath/corner/pointed/example.json'
        # fnm = '/Users/mricardo/Documents/MATLAB/optimisation/discretize/form2.json'
        fnm = '/Users/mricardo/compas_dev/me/loadpath/corner/pointed/rounded.json'
        # fnm = '/Users/mricardo/compas_dev/me/loadpath/Midsupport/topology/'+pat+'_lp.json'
        form = FormDiagram.from_json(fnm)
        # overview_forces(form)

        lp = 0
        try:
            t = form.attributes['offset']
            print('Translation of {0:.2f}'.format(t))
        except BaseException:
            t = 0.0
            print('No offset!')

        # thrust_layer = 'Thrust_grad_lp'
        # thrust_layer = str(j)+'_Fix::Thrust-'+str(i)
        thrust_layer = 'Default' #'Topology::Midsupport_'+pat+'_Thrust'
        # thrust_layer = str(j)+'_midsupport::' + str(i) + '_opt'
        # reactions_layer = 'Thrust_grad_lp'
        # points_layer = '1_Fix::Points'

        # thrust_layer = 'Thrust'
        rs.AddLayer(thrust_layer)
        # rs.AddLayer(reactions_layer)

        # i += 1

        # print(lp)

        # artist = NetworkArtist(form, layer=thrust_layer)
        # artist.clear_layer()

        #form.set_vertices_attributes(form.vertices(), {'z': 0})

        #for key, attr in form.vertices(True):
        #    form.vertex[key]['z'] = 0.0
        #    coord = form.vertex_coordinates(key)
        #    point = rs.AddPoint(coord)
        #    rs.ObjectName(point, str(form.vertex[key]['pz']))

        rs.CurrentLayer(thrust_layer)
        # rs.CurrentLayer(points_layer)

        dx = (i-1)*15
        dy = 0#-20*6 -20*j #+40 #-20 * 4

        for uv in form.edges():
            u, v = uv
            q = form.edge_attribute(uv, 'q')
            l = form.edge_length(u,v)
            if form.edge_attribute((u,v), 'is_symmetry') == False and form.edge_attribute((u,v), '_is_edge') == True and form.edge_attribute((u,v), '_is_external') == False :
                lp += q * l * l
                sp = form.vertex_coordinates(u)
                sp[0] += dx
                sp[1] += dy
                ep = form.vertex_coordinates(v)
                ep[0] += dx
                ep[1] += dy
                pz = form.vertex_attribute(u, 'pz')
                # sp[2] = form.vertex_attribute(u, 'target')
                # ep[2] = form.vertex_attribute(v, 'target')
                # sp[2] += t
                # ep[2] += t
                id = rs.AddLine(sp, ep)
                rs.ObjectName(id, str(q))
                # rs.AddTextDot(str(round(pz,2)),sp)
                # rs.ObjectColor(id, (255,255*(1002-i)/1002,0))

        rs.AddTextDot('{0:.1f}'.format(lp),[dx - 1.0, dy - 1.0, 0.0 ])
        i = i+1
        # artist = NetworkArtist(form, layer=reactions_layer)
        # artist.clear_layer()
        # rs.CurrentLayer(reactions_layer)

        # pzt=0
        # for key in form.vertices():
        #     pz = form.vertex_attribute(key,'pz')
        #     pzt+= pz
        #     if form.vertex[key]['is_fixed'] is True:
        #         node = form.vertex_coordinates(key)
        #         node[2] += t
        #         ry = form.vertex_attribute(key, '_ry')
        #         rx = form.vertex_attribute(key, '_rx')
        #         rz = form.vertex_attribute(key, '_rz', 0.0)
        #         norm = (rx ** 2 + ry ** 2 + rz ** 2) ** (1/2)
        #         print(rx,ry,rz)
        #         print(norm)
        #         if rz < 0.0 and norm > 0.0:
        #             sp = node
        #             print(sp)
        #             dz = rz/norm
        #             mult = node[2]/dz
        #             dz *= mult
        #             dx = mult* rx/norm
        #             dy = mult* ry/norm
        #             ep = [sp[0]-dx, sp[1]-dy, sp[2]-dz]
        #             id = rs.AddLine(sp, ep)
        #             rs.ObjectName(id, str(norm))
        #             rs.AddTextDot('ry: {0:.1f} / rx: {0:.1f}'.format(ry,rx),node)

        # print(pzt)

        #rs.LayerVisible('Dots', False)
        #rs.LayerVisible('Thrust', False)
