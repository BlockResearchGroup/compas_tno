
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers import view_shapes
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from numpy import array
from compas_tno.datastructures import MeshDos
from compas_plotters import MeshPlotter
import math


discretisations = [10, 20]
# degs = [10, 15, 20, 25, 30, 40]
degs = [0, 5, 10, 15, 20, 25, 30, 35, 40]
# degs = [20, ]

thk_min = 0.20
t0 = 0.50

for discretisation in discretisations:
    sols = {}

    type_formdiagram = 'cross_fd'
    type_structure = 'crossvault'
    span = span_x = span_y = 10.0
    xy_span = [[0, span_x], [0, span_y]]
    k = 1.0
    n = 4
    t = 0.0

    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0, span], [0, k*span]],
        'discretisation': discretisation,
        'fix': 'corners',
    }

    form = FormDiagram.from_library(data_diagram)

    j = 0

    z_corner_ub_strategy_middle = []
    z_corner_lb_strategy_middle = []
    z_corner_ub_strategy_ublb = []
    z_corner_lb_strategy_ublb = []
    z_corner_ub_strategy_intrados = []
    z_corner_lb_strategy_intrados = []
    z_corner_lb_analytic = []
    z_corner_ub_analytic = []

    for deg in degs:

        A = 1/math.cos(math.radians(deg))
        # print('A, deg, discretisation:', A, deg, discretisation)
        xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
        thk_min = 0.20

        ## STRATEGY ANALYTICAL:

        data_shape = {
            'type': type_structure,
            'thk': thk_min,
            'discretisation': discretisation*n,
            'xy_span': xy_span_shape,
            't': t,
        }

        vault_analytical = Shape.from_library(data_shape)
        # view_shapes(vault_analytical).show()

        ## 1. STRATEGY MIDDLE:

        xyz = array(form.vertices_attributes('xyz'))
        zt = crossvault_middle_update(xyz[:, 0], xyz[:, 1],  t,  xy_span=xy_span_shape)
        middle = MeshDos.from_mesh(form)
        i = 0
        for key in middle.vertices():
            middle.vertex_attribute(key, 'z', zt[i])
            i += 1

        vault_middle = Shape.from_middle(middle, thk=thk_min, treat_creases=True)
        # view_shapes(vault_middle).show()

        ## 2. STRATEGY OFFSET:

        t0 = 0.50

        xyz = array(form.vertices_attributes('xyz'))
        zub, zlb = crossvault_ub_lb_update(xyz[:, 0], xyz[:, 1], t0, t,  xy_span=xy_span_shape)
        intrados = MeshDos.from_mesh(form)
        extrados = MeshDos.from_mesh(form)
        middle = MeshDos.from_mesh(form)
        i = 0
        for key in middle.vertices():
            intrados.vertex_attribute(key, 'z', zlb[i])
            extrados.vertex_attribute(key, 'z', zub[i])
            middle.vertex_attribute(key, 'z', zt[i])
            i += 1

        vault_ub_lb = Shape.from_meshes(intrados, extrados, middle=middle, data={'type': 'general', 't': 0.0, 'thk': t0, 'xy_span': xy_span})
        # view_shapes(vault_ub_lb).show()

        vault_ub_lb.intrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
        vault_ub_lb.extrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
        vault_ub_lb.intrados.store_normals(correct_creases=True)
        vault_ub_lb.extrados.store_normals(correct_creases=True)

        n_offset = (t0 - thk_min)/2
        vault_ub_lb.intrados = vault_ub_lb.intrados.offset_mesh(n=n_offset, direction='up')
        vault_ub_lb.extrados = vault_ub_lb.extrados.offset_mesh(n=n_offset, direction='down')

        # view_shapes(vault_ub_lb).show()

        ## 3. STRATEGY INTRADOS OFFSET:

        zub, zlb = crossvault_ub_lb_update(xyz[:, 0], xyz[:, 1], thk_min, t,  xy_span=xy_span_shape)
        zt = crossvault_middle_update(xyz[:, 0], xyz[:, 1],  t,  xy_span=xy_span_shape)

        i = 0
        for key in form.vertices():
            form.vertex_attribute(key, 'lb', zlb[i])
            form.vertex_attribute(key, 'ub', zub[i])
            i += 1

        # form.envelope_from_shape(vault_analytical)

        intrados = MeshDos.from_mesh(form)
        i = 0
        for key in intrados.vertices():
            intrados.vertex_attribute(key, 'z', zlb[i])
            i += 1

        intrados.identify_creases_at_diagonals(xy_span=data_diagram['xy_span'])
        intrados.store_normals(correct_creases=True)
        mult = 1.0
        for key in form.fixed():
            intra_n = intrados.vertex_attribute(key, 'n')
            intrados.vertex_attribute(key, 'n', [mult*intra_n[0], mult*intra_n[1], mult*intra_n[2]])

        extrados = intrados.offset_mesh(n=thk_min, direction='up')
        middle = intrados.offset_mesh(n=thk_min/2, direction='up')

        vault_intrados = Shape.from_meshes(intrados, extrados, middle=middle, data={'type': 'general', 't': 0.0, 'thk': thk_min})
        # view_shapes(vault_intrados).show()

        diff_middle = []
        diff_ub_lb = []
        diff_intrados = []

        # view_shapes(vault_analytical).show()
        # view_shapes(vault_middle).show()
        # view_shapes(vault_ub_lb).show()
        # view_shapes(vault_intrados).show()

        print('PROBLEM')
        print('deg / discretisation:', deg, discretisation)

        print('INTRADOS')

        exit_ = True

        for key in form.vertices():
            z = form.vertex_attribute(key, 'lb')
            z_middle = vault_middle.intrados.vertex_attribute(key, 'z')[0]
            z_ub_lb = vault_ub_lb.intrados.vertex_attribute(key, 'z')[0]
            z_intrados = vault_intrados.intrados.vertex_attribute(key, 'z')[0]
            # if form.vertex_coordinates(key)[:2] == [1.0, 1.0] and exit_:
            if form.vertex_attribute(key, 'is_fixed') and exit_:
                z_corner_lb_strategy_middle.append(z_middle)
                z_corner_lb_strategy_ublb.append(z_ub_lb)
                z_corner_lb_strategy_intrados.append(z_intrados)
                z_corner_lb_analytic.append(z)
                print([form.vertex_coordinates(key)[:2], z])
                print(vault_intrados.intrados.vertex_coordinates(key))
                print(vault_ub_lb.intrados.vertex_coordinates(key))
                print(vault_middle.intrados.vertex_coordinates(key))
                print(z, z_middle, z_ub_lb, z_intrados)
                exit_ = False
            diff_middle.append(z_middle - z)
            diff_ub_lb.append(z_ub_lb - z)
            diff_intrados.append(z_intrados - z)

        print('min/max - diff middle:', min(diff_middle), max(diff_middle))
        print('min/max - diff ub_lb:', min(diff_ub_lb), max(diff_ub_lb))
        print('min/max - diff intrados:', min(diff_intrados), max(diff_intrados))

        # for mesh, diff in [[vault_middle.intrados, diff_middle], [vault_ub_lb.intrados, diff_ub_lb], [vault_intrados.intrados, diff_intrados]]:
        # for mesh, diff in [[vault_middle.intrados, diff_middle]]:
        #     plotter = MeshPlotter(mesh)
        #     plotter.draw_edges()
        #     plotter.draw_vertices(text={key: round(diff[i], 2) for i, key in enumerate(mesh.vertices())})
        #     plotter.show()

        print('-'*5)

        diff_middle = []
        diff_ub_lb = []
        diff_intrados = []

        print('EXTRADOS')

        exit_ = True

        for key in form.vertices():
            z = form.vertex_attribute(key, 'ub')
            z_middle = vault_middle.extrados.vertex_attribute(key, 'z')[0]
            z_ub_lb = vault_ub_lb.extrados.vertex_attribute(key, 'z')[0]
            z_intrados = vault_intrados.extrados.vertex_attribute(key, 'z')[0]
            # if form.vertex_coordinates(key)[:2] == [1.0, 1.0] and exit_:
            if form.vertex_attribute(key, 'is_fixed') and exit_:
                z_corner_ub_strategy_middle.append(z_middle)
                z_corner_ub_strategy_ublb.append(z_ub_lb)
                z_corner_ub_strategy_intrados.append(z_intrados)
                z_corner_ub_analytic.append(z)
                print([form.vertex_coordinates(key)[:2], z])
                print(vault_intrados.extrados.vertex_coordinates(key))
                print(vault_ub_lb.extrados.vertex_coordinates(key))
                print(vault_middle.extrados.vertex_coordinates(key))
                print(z, z_middle, z_ub_lb, z_intrados)
                print(z, z_middle, z_ub_lb, z_intrados)
                exit_ = False
            diff_middle.append(z_middle - z)
            diff_ub_lb.append(z_ub_lb - z)
            diff_intrados.append(z_intrados - z)

        print('min/max - diff middle:', min(diff_middle), max(diff_middle))
        print('min/max - diff ub_lb:', min(diff_ub_lb), max(diff_ub_lb))
        print('min/max - diff intrados:', min(diff_intrados), max(diff_intrados))

        # for mesh, diff in [[vault_middle.intrados, diff_middle], [vault_ub_lb.intrados, diff_ub_lb], [vault_intrados.intrados, diff_intrados]]:
        # for mesh, diff in [[vault_middle.intrados, diff_middle]]:
        #     plotter = MeshPlotter(mesh)
        #     plotter.draw_edges()
        #     plotter.draw_vertices(text={key: round(diff[i], 2) for i, key in enumerate(mesh.vertices())})
        #     plotter.show()

        print('-'*5)

        form.to_json('/Users/mricardo/compas_dev/me/min_thk/test/' + str(deg) + '.json')

    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(degs, z_corner_lb_analytic, label='lb_analytic')
    ax.plot(degs, z_corner_ub_analytic, label='ub_analytic')
    # ax.plot(degs, z_corner_lb_strategy_middle, label='lb_middle')
    # ax.plot(degs, z_corner_ub_strategy_middle, label='ub_middle')
    # ax.plot(degs, z_corner_lb_strategy_ublb, label='lb_ublb')
    # ax.plot(degs, z_corner_ub_strategy_ublb, label='ub_ublb')
    ax.plot(degs, z_corner_lb_strategy_intrados, label='lb_intrados')
    ax.plot(degs, z_corner_ub_strategy_intrados, label='ub_intrados')
    ax.legend()
    plt.show()
