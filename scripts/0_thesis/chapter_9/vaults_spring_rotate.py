from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas.geometry import norm_vector
from compas.geometry import Line, Point, Vector
from compas_tno.utilities.form import slide_pattern_inwards
from compas.colors import Color
import math
from numpy import array

discretisation = 14
spr_angle = 30.0
R_over_L = 0.5

L_form = 10.0
xf = L_form
x0 = 0.0
xc = yc = (x0 + xf)/2
printout = False
starting = 'loadpath'
solver = 'IPOPT'

cos = math.cos(math.radians(spr_angle))
fact = 2 * (R_over_L * (cos - 1) + 1/2)
L = L_form/fact
R = L * R_over_L
Ldiff = L - L_form
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

results = {}

vault = Shape.create_crossvault(discretisation=discretisation*2, xy_span=xyspan_shape)
vault.ro = 0.1

ro = 45.0
delta = 0.714

for lambd in [0.0]:

    form = FormDiagram.create_cross_form(discretisation=discretisation)
    # form = FormDiagram.create_fan_form(discretisation=discretisation)
    # form = FormDiagram.create_parametric_form(lambd=lambd, discretisation=discretisation)

    # form = FormDiagram.create_cross_with_diagonal(discretisation=discretisation)
    # slide_pattern_inwards(form, delta=delta)
    # form_base = FormDiagram.create_cross_form(discretisation=discretisation)
    # slide_pattern_inwards(form_base, delta=delta)

    starting = 'loadpath'

    print('-'*10, ' Analysis for lambda= ', lambd)

    results[lambd] = {}

    # for i in range(36):
    for i in [19]:

        phi = 10 * i
        results[lambd][phi] = {}

        print('-'*10, ' Analysis for i/phi= ',i, phi)

        key_sup = None
        vector_supports = []
        lines = []
        plots_vectors = []
        for key in form.vertices_where({'is_fixed': True}):
            x, y, z = form.vertex_coordinates(key)
            dXbi = [0, 0, 0]
            if x > xc and y > yc:
                key_sup = key
                dXbi = [math.cos(math.radians(ro)) * math.cos(math.radians(-phi)), math.cos(math.radians(ro)) * math.cos(math.radians(-phi)),  math.sin(math.radians(-phi))]
                print('Norm of vector:', dXbi, norm_vector(dXbi))
                start = [x, y, z]
                end = [x + dXbi[0], y + dXbi[1], z + dXbi[2]]
                lines.append(Line(start, end))
                plots_vectors.append([Point(*start), Vector(*dXbi)])
            vector_supports.append(dXbi)

        dXb = array(vector_supports)
        # print(dXb)

        problem: Analysis = Analysis.create_compl_energy_analysis(form, vault,
                                                                  printout=printout,
                                                                  starting_point=starting,
                                                                  support_displacement=dXb,
                                                                  solver=solver)
        problem.apply_selfweight()
        # problem.apply_selfweight_from_pattern(form_base)
        problem.apply_envelope()
        problem.set_up_optimiser()
        problem.run()

        swt = abs(form.lumped_swt())
        # print('SWT:', swt)

        fopt = problem.optimiser.fopt
        rx = form.vertex_attribute(key_sup, '_rx')
        ry = form.vertex_attribute(key_sup, '_ry')
        V = abs(form.vertex_attribute(key_sup, '_rz'))
        T = math.sqrt(rx**2 + ry**2)

        results[lambd][phi]['fopt'] = ' '
        results[lambd][phi]['V/W'] = ' '
        results[lambd][phi]['T/W'] = ' '

        if problem.optimiser.exitflag == 0:
            print('Solved for i=', i)
            starting = 'current'
            results[lambd][phi]['fopt'] = fopt
            results[lambd][phi]['V/W'] = V/swt
            results[lambd][phi]['T/W'] = T/swt

        else:
            print('No sol for i=', i)
            starting = 'loadpath'

        view: Viewer = Viewer(form)
        # view.scale_edge_thickness(10.0)
        # view.draw_thrust()
        view.draw_shape()
        # view.draw_cracks()
        # for base, vector in plots_vectors:
        #     base_sup = Point(*form.vertex_coordinates(key_sup))
        #     view.draw_vector(vector, base_sup, color=Color.black())
        view.show()

    print(results)

for lambd in results:
    print('-'*20)
    print('lambd:', lambd)
    for phi in results[lambd]:
        print(phi, results[lambd][phi]['fopt'], results[lambd][phi]['V/W'], results[lambd][phi]['T/W'])

