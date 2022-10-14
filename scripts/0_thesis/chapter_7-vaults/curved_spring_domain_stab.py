from ast import For
from compas_tno import analysis
from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities.form import slide_pattern_inwards
from compas_tno.utilities.loads import apply_selfweight_from_shape
import math
import os

form = FormDiagram.create_cross_form()
vault = Shape.create_pointedcrossvault()

R_over_L = 0.2
# lambd = 0.8
thk_0 = 0.5

# tests

printout = False

results = {}

for spr_angle in [20]:

    results[spr_angle] = {}

    for delta in [0.3, 0.4, 0.5, 0.6, 0.7]:  # 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0

        print('-'*10, ' Analysis for delta= ', delta)

        results[spr_angle][delta] = {}

        discr = 14
        L_form = 10.0
        xf = L_form
        x0 = 0.0
        xc = yc = (x0 + xf)/2
        xyspan = [[x0, xf], [x0, xf]]

        form = FormDiagram.create_cross_with_diagonal(xy_span=xyspan, discretisation=discr)
        slide_pattern_inwards(form, delta=delta)

        form_base = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=discr)
        slide_pattern_inwards(form_base, delta=delta)

        i = 0

        for R_over_L in [0.79422]:  # 0.61147, 0.708, 0.79422

            print('-'*10, ' Analysis for R_over_L= ', R_over_L)

            folder = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/curved_fd/delta_{}/spr_{}/stabiltydomain/'.format(delta, spr_angle)
            os.makedirs(folder, exist_ok=True)
            # folder_minthk = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/curved_fd/delta_{}/spr_{}/'.format(delta, spr_angle)
            # min_thk = 'curved_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}.json'.format(discr, delta, spr_angle, R_over_L)
            # form = FormDiagram.from_json(folder_minthk + min_thk)

            results[spr_angle][delta][R_over_L] = {}

            cos = math.cos(math.radians(spr_angle))
            fact = 2 * (R_over_L * (cos - 1) + 1/2)
            L = L_form/fact
            R = L * R_over_L
            Ldiff = L - L_form
            xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

            hc = math.sqrt(R**2 - (R - L/2)**2)

            if R_over_L == 0.61147:
                thks = [0.15]  # [::-1]  0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20
            if R_over_L == 0.708:
                thks = [0.25]  # [::-1]
            if R_over_L == 0.79422:
                thks = [0.30]  # [::-1]

            for thk in thks:

                if i == 0:
                    starting = 'loadpath'
                else:
                    starting = 'current'

                results[spr_angle][delta][R_over_L][str(round(thk, 2))] = {}
                results[spr_angle][delta][R_over_L][str(round(thk, 2))]['min'] = ' '
                results[spr_angle][delta][R_over_L][str(round(thk, 2))]['max'] = ' '

                vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk)
                vault.ro = 0.1  # really important to IPOPT solution

                apply_selfweight_from_shape(form_base, vault)

                problem: Analysis = Analysis.create_minthrust_analysis(form, vault, printout=printout, solver='SLSQP', starting_point=starting)
                # problem.apply_selfweight()
                problem.apply_selfweight_from_pattern(form_base)
                problem.apply_envelope()
                problem.set_up_optimiser()
                problem.run()

                i += 1

                T_over_W = abs(form.thrust()/form.lumped_swt())

                if problem.optimiser.exitflag == 0:
                    results[spr_angle][delta][R_over_L][str(round(thk, 2))]['min'] = T_over_W
                    title = 'curved_fd_discr_{}_minthk_delta_{}_spr_{}_R_over_L_{}_thk_{}_min.json'.format(discr, delta, spr_angle, R_over_L, thk)
                    path = folder+title
                    form.to_json(path)
                    print('Saved:', path)
                else:
                    i = 0
                    break

            for thk in thks:

                try:  # load min thk
                    title = 'curved_fd_discr_{}_minthk_delta_{}_spr_{}_R_over_L_{}_thk_{}_min.json'.format(discr, delta, spr_angle, R_over_L, thk)
                    form = FormDiagram.from_json(folder+title)
                    print('Trying to load:', folder+title)
                    starting = 'current'
                except:
                    print('Did not load:', folder+title)
                    if i == 0:
                        starting = 'loadpath'
                    else:
                        starting = 'current'

                # results[spr_angle][delta][R_over_L][str(round(thk, 2))] = {}  # remove
                # results[spr_angle][delta][R_over_L][str(round(thk, 2))]['min'] = ' '  # remove
                results[spr_angle][delta][R_over_L][str(round(thk, 2))]['max'] = ' '

                vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk)
                vault.ro = 0.1  # really important to IPOPT solution

                apply_selfweight_from_shape(form_base, vault)

                problem: Analysis = Analysis.create_maxthrust_analysis(form, vault, printout=printout, solver='SLSQP', starting_point=starting)
                # problem.apply_selfweight()  # check!!!
                problem.apply_selfweight_from_pattern(form_base)
                problem.apply_envelope()
                problem.set_up_optimiser()
                problem.run()
                i += 1

                T_over_W = abs(form.thrust()/form.lumped_swt())

                if problem.optimiser.exitflag == 0:
                    results[spr_angle][delta][R_over_L][str(round(thk, 2))]['max'] = T_over_W
                    title = 'curved_fd_discr_{}_minthk_delta_{}_spr_{}_R_over_L_{}_thk_{}_max.json'.format(discr, delta, spr_angle, R_over_L, thk)
                    path = folder+title
                    form.to_json(path)
                    print('Saved:', path)
                else:
                    i = 0
                    # break

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
            thks = results[spr_angle][lambd][R_over_L]
            for thk in thks:
                try:
                    print(thk, results[spr_angle][lambd][R_over_L][thk]['min'], results[spr_angle][lambd][R_over_L][thk]['max'])
                except:
                    print(thk, results[spr_angle][lambd][R_over_L][thk]['min'])

