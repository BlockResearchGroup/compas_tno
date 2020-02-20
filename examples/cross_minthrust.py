from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form

from compas_tno.utilities.constraints import set_cross_vault_heights
from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters.plotters import plot_form

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'cross_fd'
    objective = 'max'
    thck = 0.3044
    decrease = 0.0001
    exitflag = 0

    while exitflag == 0:

        # Create Vault from one of the patterns Fan/Grid with the dimensions

        print('/ ----------- THCK: ', thck)

        x_span = 7.5           # Dont forget to ajust mesures accordingly
        y_span = 10.0
        example = 'rectangular/7,5x10/'

        if type_fd == 'cross_fd':
            divisions = 20
            form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
        if type_fd == 'fan_fd':
            divisions = 16
            form = create_fan_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
        if type_fd == 'mixed_fd':
            divisions = 16
            PATH = '/Users/mricardo/compas_dev/me/minmax/cross/' + example + type_fd + '/' + type_fd + '_discr_'+ str(divisions)
            file_initial = PATH + '_lp.json'
            form = FormDiagram.from_json(file_initial)

        PATH = '/Users/mricardo/compas_dev/me/minmax/cross/' + example + type_fd + '/' + type_fd + '_discr_'+ str(divisions)
        thck_save = thck + 0.0001
        file_save = PATH + '_' + objective + '_t=' + str(int(thck_save*1000)) + '.json'
        print('File to be saved: ',file_save)
        # file_initial = PATH + '_lp.json'
        thck_initial = thck + decrease + 0.0001
        # thck_initial = 0.420001
        objective_initial = 'min' # objective
        file_initial = PATH + '_' + objective_initial + '_t=' + str(int((thck_initial*1000))) + '.json'
        print('File initial: ',file_initial)
        form = FormDiagram.from_json(file_initial)

        # Set Constraints for Cross_Vaults

        tol = 0.0
        # thck = 0.31
        form = set_cross_vault_heights(form, xy_span = [[0.0 - tol,x_span + tol],[0.0 - tol,y_span + tol]], thk = thck, b = 5.0, set_heights=False, ub_lb = True, update_loads = True)

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
        indset = form.attributes['indset']

        # Maximum or Minimum Thrusts

        solver = 'pyOpt-SLSQP'
        # solver = 'slsqp'

        fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                            printout=print_opt,
                                            find_inds=True,
                                            indset=indset,
                                            tol=0.001,
                                            translation = translation,
                                            objective=objective,
                                            bmax = False,
                                            summary=True)

        if exitflag == 0:
            print('File saved to: ',file_save)
            form.to_json(file_save)
            overview_forces(form)
            reactions(form)
            # plot_form(form, show_q = False).show()

        for key in form.vertices_where({'is_fixed': True}):
            rx = round(form.get_vertex_attribute(key, 'rx'),3)
            ry = round(form.get_vertex_attribute(key, 'ry'),3)
            zb = round(form.get_vertex_attribute(key,'z'),3)
            print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
            break
        qmax = round(max(qopt).item(),3)
        fopt = round(fopt,3)
        print('zb: {0:.3f}'.format(zb))
        print('fopt: {0:.3f}'.format(fopt))
        print('qmax: {0:.3f}'.format(max(qopt).item()))

        thck = round(thck - decrease,4)


