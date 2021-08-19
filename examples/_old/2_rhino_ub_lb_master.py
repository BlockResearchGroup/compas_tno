from compas_tna.diagrams import FormDiagram
from compas_rhino.artists import NetworkArtist
from compas.utilities import geometric_key
from compas.geometry import distance_point_point
from compas.utilities import i_to_rgb
# from compas_tno.diagrams.form import overview_forces
# from compas.utilities import XFunc

import rhinoscriptsyntax as rs
import json
import math

# Mesh

for t in [0.14]:
    for objective in ['min']:

        fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/fan_fd/fan_fd_discr_16_'+ objective + '_t='+ str(int(t*100)) +'.json'

        fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular-rollers/7,5x10/cross_fd/cross_fd_discr_20_rol-all_min_t=14.json'

        # fnm = '/Users/mricardo/compas_dev/me/minmax/dome/flower/flower_discr_8_16_'+ objective + '_t='+ str(int(t*100)) +'.json'
        print('Load: ', fnm)
        form = FormDiagram.from_json(fnm)
        k_i = form.key_index()
        i_k = form.index_key()

        radius_circlus = 0.10

        master = 'cross_rol_t=-' + str(t) + '_' + objective

        try:
            cracks_lb, cracks_ub = form.attributes['cracks']
        except BaseException:
            cracks_lb = []
            cracks_ub = []
        fs = []
        radius_max = 0.30
        # Usually 0.50 for Cross Form-Diagram

        thrust_layer = master + '::Thrust'
        rs.AddLayer(thrust_layer)
        rs.CurrentLayer(thrust_layer)
        artist = NetworkArtist(form, layer=thrust_layer)
        artist.clear_layer()
        lp = 0.0
        for u, v in form.edges():
            q = form.edge_attribute((u,v), 'q')
            l = form.edge_length(u,v)
            fs.append(q*l)
            sp = form.vertex_coordinates(u)
            ep = form.vertex_coordinates(v)
            id = rs.AddLine(sp, ep)
            rs.ObjectName(id, str(q))
            lp += q * l * l
        rs.AddTextDot('{0:.1f}'.format(lp),[- 1.0, - 1.0, 0.0 ])

        thrust_limits = master + '::Limits'
        rs.AddLayer(thrust_limits)
        rs.CurrentLayer(thrust_limits)
        artist = NetworkArtist(form, layer=thrust_limits)
        artist.clear_layer()
        lp = 0.0
        tol_crack = 0.01
        for key in form.vertices():
            x, y, z = form.vertex_coordinates(key)
            try:
                lb = form.vertex_attribute(key, 'lb')
                if abs(z - lb) < tol_crack:
                    id = rs.AddPoint([x,y,z])
                    rs.ObjectColor(id, color = (0,0,200))
                    id = rs.AddSphere([x,y,z], radius_circlus)
                    rs.ObjectColor(id, color = (0,0,200))
            except BaseException:
                pass
            try:
                ub = form.vertex_attribute(key, 'ub')
                if abs(z - ub) < tol_crack:
                    id = rs.AddPoint([x,y,z])
                    rs.ObjectColor(id, color = (200,0,0))
                    id = rs.AddSphere([x,y,z], radius_circlus)
                    rs.ObjectColor(id, color = (200,0,0))
            except BaseException:
                pass
        if cracks_lb:
            for i in cracks_lb:
                rs.AddTextDot('lb',form.vertex_coordinates(k_i[i]))
        if cracks_ub:
            for i in cracks_ub:
                rs.AddTextDot('ub',form.vertex_coordinates(k_i[i]))



        pipes_color = master + '::Pipes-Color'
        rs.AddLayer(pipes_color)
        rs.CurrentLayer(pipes_color)
        artist = NetworkArtist(form, layer=pipes_color)
        artist.clear_layer()
        rs.CurrentLayer(pipes_color)
        for u, v in form.edges_where({'_is_edge': True}):
            l = form.edge_length(u, v)
            q = form.edge_attribute((u, v), 'q')
            sp = form.vertex_coordinates(u)
            ep = form.vertex_coordinates(v)
            if q > 0.001:
                id = rs.AddLine(sp, ep)
                f = math.fabs(q * l)
                coef = f/max(fs)
                r, g, b = i_to_rgb(coef)
                # print(coef)
                pipe = rs.AddPipe(id, 0, 0.06)
                rs.ObjectColor(pipe, color=(r,g,b))

        rs.CurrentLayer('Default')
        rs.LayerVisible(thrust_layer, False)
        rs.LayerVisible(thrust_limits, False)
        rs.LayerVisible(pipes_color, False)
        print(max(fs))
