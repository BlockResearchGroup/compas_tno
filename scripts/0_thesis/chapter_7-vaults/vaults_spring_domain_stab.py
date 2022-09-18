from compas_tno import analysis
from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
import math
import os

form = FormDiagram.create_cross_form()
vault = Shape.create_pointedcrossvault()

R_over_L = 0.2
lambd = 0.8
thk_0 = 0.5

results = {}

for spr_angle in [20]:

    results[spr_angle] = {}

    for lambd in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:  # 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # for lambd in [0.9]:

        print('-'*10, ' Analysis for lambda= ', lambd)

        results[spr_angle][lambd] = {}

        discr = 14
        L_form = 10.0
        xf = L_form
        x0 = 0.0
        xc = yc = (x0 + xf)/2
        xyspan = [[x0, xf], [x0, xf]]

        form = FormDiagram.create_parametric_form(xy_span=xyspan, lambd=lambd, discretisation=discr)

        i = 0

        for R_over_L in [0.708]:  # 0.61147, 0.708, 0.79422

            print('-'*10, ' Analysis for R_over_L= ', R_over_L)

            results[spr_angle][lambd][R_over_L] = {}

            cos = math.cos(math.radians(spr_angle))
            fact = 2 * (R_over_L * (cos - 1) + 1/2)
            L = L_form/fact
            R = L * R_over_L
            Ldiff = L - L_form
            xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

            hc = math.sqrt(R**2 - (R - L/2)**2)

            # if R_over_L == 0.61147:
            #     thks = [0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]
            if R_over_L == 0.708:
                thks = [0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20]
            if R_over_L == 0.708:
                thks = [0.25]
            # if R_over_L == 0.79422:
            #     thks = [0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20]

            for thk in thks:

                if i == 0:
                    starting = 'loadpath'
                else:
                    starting = 'current'

                results[spr_angle][lambd][R_over_L][str(round(thk, 2))] = {}
                results[spr_angle][lambd][R_over_L][str(round(thk, 2))]['min'] = ' '

                vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk)
                # vault.ro = 0.1  # really important to IPOPT solution

                problem: Analysis = Analysis.create_minthrust_analysis(form, vault, printout=False, solver='SLSQP', starting_point=starting)
                problem.apply_selfweight()
                problem.apply_envelope()
                problem.set_up_optimiser()
                problem.run()

                i += 1

                T_over_W = abs(form.thrust()/form.lumped_swt())

                folder = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/parametric_fd/lambd_{}/spr_{}/stabiltydomain/'.format(lambd, spr_angle)
                os.makedirs(folder, exist_ok=True)

                if problem.optimiser.exitflag == 0:
                    results[spr_angle][lambd][R_over_L][str(round(thk, 2))]['min'] = T_over_W
                    title = 'parametric_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}_thk_{}_min.json'.format(discr, lambd, spr_angle, R_over_L, thk)
                    path = folder+title
                    form.to_json(path)
                    print('Saved:', path)
                else:
                    i = 0
                    break

            for thk in thks:

                if i == 0:
                    starting = 'loadpath'
                else:
                    starting = 'current'

                results[spr_angle][lambd][R_over_L][str(round(thk, 2))]['max'] = ' '

                vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk)
                # vault.ro = 0.1  # really important to IPOPT solution

                problem: Analysis = Analysis.create_maxthrust_analysis(form, vault, printout=False, solver='SLSQP', starting_point=starting)
                problem.apply_selfweight()
                problem.apply_envelope()
                problem.set_up_optimiser()
                problem.run()
                i += 1

                T_over_W = abs(form.thrust()/form.lumped_swt())

                if problem.optimiser.exitflag == 0:
                    results[spr_angle][lambd][R_over_L][str(round(thk, 2))]['max'] = T_over_W
                    title = 'parametric_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}_thk_{}_max.json'.format(discr, lambd, spr_angle, R_over_L, thk)
                    path = folder+title
                    form.to_json(path)
                    print('Saved:', path)
                else:
                    i = 0
                    break

            print(results)

        print(results)

print(results)

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
            print('>'*10)
            print(R_over_L)
            for thk in results[spr_angle][lambd][R_over_L]:
                print(thk, results[spr_angle][lambd][R_over_L][thk]['min'], results[spr_angle][lambd][R_over_L][thk]['max'])

