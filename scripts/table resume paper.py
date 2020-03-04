
from compas_tna.diagrams import FormDiagram
from compas_tno.algorithms.equilibrium import reactions
from compas_tno.plotters import plot_form
import math

# Resume data for Rectangular Cross-Vault.

# Cross Diagram

print('\n\n/----------- CROSS DIAGRAM -----------')

for t in [0.50, 0.305]:
    for objective in ['min', 'max']:
        try:
            fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/cross_fd/cross_fd_discr_20_'+ objective + '_t='+ str(int(t*1000)) +'.json'
            form = FormDiagram.from_json(fnm)
        except:
            fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/cross_fd/cross_fd_discr_20_'+ objective + '_t='+ str(int(t*100)) +'.json'
            form = FormDiagram.from_json(fnm)
        print('\n\n/------------- Load: ', t, objective, '\n', fnm)
        form = FormDiagram.from_json(fnm)
        q = [form.edge_attribute((u,v), 'q') for u,v in form.edges()]
        f = [form.edge_attribute((u,v), 'q')*form.edge_length(u,v) for u,v in form.edges()]
        zb = [form.vertex_attribute(key, 'z') for key in form.vertices_where({'is_fixed': True})]
        print('Q max/min:', round(max(q),3), round(min(q),3))
        print('f max/min:', round(max(f),3), round(min(f),3))
        print('zb max/min:', round(max(zb),3), round(min(zb),3))
        try:
            print('Objective: ', round(form.attributes['fopt'],1))
            print('Type-opt: ', form.attributes['objective'])
        except:
            reactions(form, plot= False)
            fopt = 0
            for key in form.vertices_where({'is_fixed': True}):
                rx = form.vertex_attribute(key, 'rx')
                ry = form.vertex_attribute(key, 'ry')
                R = math.sqrt(rx**2 + ry**2)
                fopt += R
            form.attributes['fopt'] = fopt
            print('Objective: ', round(form.attributes['fopt'],1))
            print('Type-opt: ', objective)

# Fan Diagram

print('\n\n/----------- FAN DIAGRAM -----------')

for t in [0.50, 0.418]:
    for objective in ['min', 'max']:
        try:
            fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/fan_fd/fan_fd_discr_16_'+ objective + '_t='+ str(int(t*1000)) +'.json'
            form = FormDiagram.from_json(fnm)
        except:
            fnm = '/Users/mricardo/compas_dev/me/minmax/cross/rectangular/7,5x10/fan_fd/fan_fd_discr_16_'+ objective + '_t='+ str(int(t*100)) +'.json'
            form = FormDiagram.from_json(fnm)
        print('\n\n/------------- Load: ', t, objective, '\n', fnm)
        q = [form.edge_attribute((u,v), 'q') for u,v in form.edges()]
        f = [form.edge_attribute((u,v), 'q')*form.edge_length(u,v) for u,v in form.edges()]
        zb = [form.vertex_attribute(key, 'z') for key in form.vertices_where({'is_fixed': True})]
        print('Q max/min:', round(max(q),3), round(min(q),3))
        print('f max/min:', round(max(f),3), round(min(f),3))
        print('zb max/min:', round(max(zb),3), round(min(zb),3))
        print('Objective: ', round(form.attributes['fopt'],1))
        print('Type-opt: ', form.attributes['objective'])

# Resume data for Dome.

# Radial Form Diagram

print('\n\n/----------- RADIAL DIAGRAM -----------')

for t in [0.30, 0.08]:
    for objective in ['min', 'max']:
        fnm = '/Users/mricardo/compas_dev/me/minmax/dome/r=5/radial_discr_8_16_'+ objective + '_t='+ str(int(t*100)) +'.json'
        print('\n\n/------------- Load: ', t, objective, '\n', fnm)
        form = FormDiagram.from_json(fnm)
        q = [form.edge_attribute((u,v), 'q') for u,v in form.edges()]
        f = [form.edge_attribute((u,v), 'q')*form.edge_length(u,v) for u,v in form.edges()]
        zb = [form.vertex_attribute(key, 'z') for key in form.vertices_where({'is_fixed': True})]
        print('Q max/min:', round(max(q),3), round(min(q),3))
        print('f max/min:', round(max(f),3), round(min(f),3))
        print('zb max/min:', round(max(zb),3), round(min(zb),3))
        print('Objective: ', round(form.attributes['fopt'],1))
        print('Type-opt: ', form.attributes['objective'])


# Flower Form Diagram

print('\n\n/----------- FLOWER DIAGRAM -----------')

for t in [0.30, 0.25]:
    for objective in ['min', 'max']:
        fnm = '/Users/mricardo/compas_dev/me/minmax/dome/flower/flower_discr_8_16_'+ objective + '_t='+ str(int(t*100)) +'.json'
        print('\n\n/------------- Load: ', t, objective, '\n', fnm)
        form = FormDiagram.from_json(fnm)
        q = [form.edge_attribute((u,v), 'q') for u,v in form.edges()]
        f = [form.edge_attribute((u,v), 'q')*form.edge_length(u,v) for u,v in form.edges()]
        zb = [form.vertex_attribute(key, 'z') for key in form.vertices_where({'is_fixed': True})]
        print('Q max/min:', round(max(q),3), round(min(q),3))
        print('f max/min:', round(max(f),3), round(min(f),3))
        print('zb max/min:', round(max(zb),3), round(min(zb),3))
        print('Objective: ', round(form.attributes['fopt'],1))
        print('Type-opt: ', form.attributes['objective'])
