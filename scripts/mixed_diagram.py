from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form

from compas_tno.utilities.constraints import set_cross_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters import plot_form

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'mixed_fd'
    objective = 'min'
    thck = 0.50

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    x_span = 10.0
    y_span = 10.0
    example = 'square/'

    divisions = 16
    PATH = '/Users/mricardo/compas_dev/me/minmax/cross/' + example + type_fd + '/' + type_fd + '_discr_'+ str(divisions)
    file_initial = PATH + '_lp.json'
    file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '.json'

    if type_fd == 'cross_fd':
        divisions = 16
        form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
    elif type_fd == 'fan_fd':
        divisions = 16
        form = create_fan_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
    else:
        file_begin = PATH + '.json'
        form = FormDiagram.from_json(file_begin)
        # plot_form(form, show_q = False).show()

    # Set Constraints for Cross_Vaults

    # form = set_cross_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], thk = thck, b = 5.0, set_heights=False, ub_lb = True, update_loads = True)

    # Initial parameters

    translation = True
    bounds_width = 5.0
    use_bounds = False
    qmax = 100
    qmin = -1e-6
    indset = None
    print_opt = True

    # Convex Optimisation to find good starting point. Save the starting point, and can load it later if wanted

    # fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
    #                                         printout=print_opt,
    #                                         find_inds=True,
    #                                         tol=0.01,
    #                                         objective='loadpath',
    #                                         indset=indset)

    # form.to_json(file_initial)

    # print(sdud)

    form = FormDiagram.from_json(file_initial)
    indset = form.attributes['indset']

    # Maximum or Minimum Thrusts

    solver = 'slsqp'

    fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                        printout=print_opt,
                                        find_inds=True,
                                        indset=indset,
                                        tol=0.01,
                                        translation = translation,
                                        objective=objective,
                                        bmax = True,
                                        summary=print_opt)

    print('File saved to: ',file_save)
    form.to_json(file_save)
    overview_forces(form)
    # plot_form(form, show_q = False).show()
    reactions(form)

    for key in form.vertices_where({'is_fixed': True}):
        rx = round(form.vertex_attribute(key, 'rx'),3)
        ry = round(form.vertex_attribute(key, 'ry'),3)
        zb = round(form.vertex_attribute(key,'z'),3)
        print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
        break
    qmax = round(max(qopt).item(),3)
    fopt = round(fopt,3)
    print('zb: {0:.3f}'.format(zb))
    print('fopt: {0:.3f}'.format(fopt))
    print('qmax: {0:.3f}'.format(max(qopt).item()))
