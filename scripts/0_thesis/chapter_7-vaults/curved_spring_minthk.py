from compas_tno import analysis
from compas_tno.analysis.analysis import Analysis
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.utilities.loads import apply_selfweight_from_shape
from compas_tno.viewers import Viewer
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities.form import slide_pattern_inwards
import math
import csv
import os

form = FormDiagram.create_cross_form()
vault = Shape.create_pointedcrossvault()

R_over_L = 0.2
# lambd = 0.8
thk_0 = 0.5
printout = True

results = {}

# for lambd in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
for delta in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    # for lambd in [0.9]:

    print('-'*10, ' Analysis for delta= ', delta)

    results[delta] = {}

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
    for spr_angle in [20]:
    # for spr_angle in [40]:

        results[delta][spr_angle] = {}

        # 0.70, 0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795
        for R_over_L in [0.61147, 0.708, 0.79422]:  # 0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0

            results[delta][spr_angle][R_over_L] = {}

            starting = 'loadpath'

            cos = math.cos(math.radians(spr_angle))
            fact = 2 * (R_over_L * (cos - 1) + 1/2)
            L = L_form/fact
            R = L * R_over_L
            Ldiff = L - L_form
            xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

            hc = math.sqrt(R**2 - (R - L/2)**2)

            vault = Shape.create_pointedcrossvault(xy_span=xyspan_shape, discretisation=discr*2, hc=hc, thk=thk_0)
            vault.ro = 0.1  # really important to IPOPT solution

            apply_selfweight_from_shape(form_base, vault)

            # pz = {}
            # for key in form_base.vertices():
            #     pz[key] = round(form_base.vertex_attribute(key, 'pz')*10, 2)
            # plot = TNOPlotter(form_base)
            # plot.draw_form(scale_width=False)
            # plot.draw_vertexlabels(text=pz)
            # plot.show()

            problem: Analysis = Analysis.create_minthk_analysis(form, vault, printout=printout, solver='SLSQP', starting_point=starting)
            # problem.apply_selfweight()
            problem.apply_selfweight_from_pattern(form_base)
            problem.apply_envelope()

            # pz = {}
            # for key in form.vertices():
            #     pz[key] = round(form.vertex_attribute(key, 'pz')*10, 2)
            # plot = TNOPlotter(form)
            # plot.draw_form(scale_width=False)
            # plot.draw_vertexlabels(text=pz)
            # plot.show()

            problem.set_up_optimiser()
            problem.run()

            T_over_W = abs(form.thrust()/form.lumped_swt())

            i += 1

            minthk = problem.optimiser.fopt
            t_over_s = minthk/L_form

            folder = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/curved_fd/delta_{}/spr_{}'.format(delta, spr_angle)
            os.makedirs(folder, exist_ok=True)
            title = 'curved_fd_discr_{}_minthk_lambd_{}_spr_{}_R_over_L_{}.json'.format(discr, delta, spr_angle, R_over_L)
            save_form = os.path.join(folder, title)

            if problem.optimiser.exitflag == 0:
                results[delta][spr_angle][R_over_L]['t_over_s'] = t_over_s
                results[delta][spr_angle][R_over_L]['T_over_W'] = T_over_W
                form.to_json(save_form)
                print('Saved to:', save_form)
            else:
                results[delta][spr_angle][R_over_L]['t_over_s'] = ' '
                results[delta][spr_angle][R_over_L]['T_over_W'] = ' '

            # shape_nice = Shape.from_formdiagram_and_attributes(form)

            # view: Viewer = Viewer(form, show_grid=False)
            # view.scale_edge_thickness(10.0)
            # view.draw_thrust()
            # view.draw_shape()
            # view.draw_cracks()
            # view.show()

    print(results)

print(results)

for lambd in results:

    print('-'*20)
    print('-'*20)
    print('delta=', lambd)

    for spr_angle in results[lambd]:
        print('-'*20)
        print('spr_angle=', spr_angle)

        for R_over_L in results[lambd][spr_angle]:
            print(R_over_L, results[lambd][spr_angle][R_over_L]['t_over_s'], results[lambd][spr_angle][R_over_L]['T_over_W'])
