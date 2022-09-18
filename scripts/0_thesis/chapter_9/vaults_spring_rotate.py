from compas_tno import analysis
from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
import math

form = FormDiagram.create_cross_form()
vault = Shape.create_pointedcrossvault()

R_over_L = 0.2
lambd = 0.8
thk_0 = 0.5
spr_angle = 30

phi = {}

for spr_angle in [30]:

    results[spr_angle] = {}

    for lambd in [0.0]:
        # for lambd in [0.9]:

        print('-'*10, ' Analysis for lambda= ', lambd)

        results[spr_angle][lambd] = {}

        discr = 14
        L_form = 10.0
        xf = L_form
        x0 = 0.0
        xc = yc = (x0 + xf)/2
        xyspan = [[x0, xf], [x0, xf]]

        # form = FormDiagram.create_parametric_form(xy_span=xyspan, lambd=lambd, discretisation=discr)
        form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=discr)

        i = 0
        for R_over_L in [0.5]:
        # for R_over_L in [0.5]:

            if i == 0:
                starting = 'loadpath'
            else:
                starting = 'current'

            alpha = 1/math.cos(math.radians(spr_angle))
            L = xf * alpha
            R = L * R_over_L
            Ldiff = L - xf
            xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

            hc = math.sqrt(R**2 - (R - L/2)**2)

            vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk_0)
            vault.ro = 0.1  # really important to IPOPT solution

            vector_supports = []
            lines = []
            plots_vectors = []
            for key in form.vertices_where({'is_fixed': True}):
                x, y, z = form.vertex_coordinates(key)
                dXbi = [0, 0, 0]
                if x > xc and y < yc:
                    dXbi = normalize_vector([1, -1, 0])
                if norm_vector(dXbi) > 1e-3:
                    start = [x, y, z]
                    end = [x + dXbi[0], y + dXbi[1], z + dXbi[2]]
                    lines.append(Line(start, end))
                    plots_vectors.append([Point(*start), Vector(*dXbi)])
                vector_supports.append(dXbi)



            problem: Analysis = Analysis.create_minthk_analysis(form, vault, printout=False, solver='IPOPT', starting_point=starting)
            problem.apply_selfweight()
            problem.apply_envelope()
            problem.set_up_optimiser()
            problem.run()

            i += 1

            minthk = problem.optimiser.fopt
            t_over_s = minthk/L_form

            if problem.optimiser.exitflag == 0:
                results[spr_angle][lambd][R_over_L] = t_over_s

# shape_nice = Shape.from_formdiagram_and_attributes(form)

# view: Viewer = Viewer(form, show_grid=False)
# view.draw_thrust()
# view.draw_shape()
# view.draw_cracks()
# view.show()

for spr_angle in results:
    print('-'*20)
    print('-'*20)
    print(spr_angle)
    for lambd in results[spr_angle]:
        print('-'*20)
        print(lambd)
        for R_over_L in results[spr_angle][lambd]:
            print(R_over_L, results[spr_angle][lambd][R_over_L])

