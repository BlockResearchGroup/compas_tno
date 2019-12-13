from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form

from compas_tno.utilities.constraints import set_cross_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_tno.plotters.plotters import plot_form

from compas_tno.utilities.constraints import create_cracks
from compas_tno.utilities.constraints import rollers_on_openings


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    type_fd = 'cross_fd'
    objective = 'min_cracks'
    rollers = False

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    x_span = 10.0
    y_span = 10.0

    if type_fd == 'cross_fd':
        divisions = 24
        xy_span = [[0.0,x_span],[0.0,y_span]]
        form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
    if type_fd == 'fan_fd':
        divisions = 16
        xy_span = [[0.0,x_span],[0.0,y_span]]
        form = create_fan_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)

    PATH = '/Users/mricardo/compas_dev/me/minmax/cross/square/'+ type_fd + '/' + type_fd + '_discr_'+ str(divisions)
    file_initial = PATH + '_lp.json'
    file_save = PATH + '_' + objective + '.json'

    # Set Constraints for Cross_Vaults

    form = set_cross_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], thk = 0.5, b = 5.0, set_heights=False, ub_lb = True, update_loads = True)

    # Initial parameters

    translation = form.attributes['tmax']
    qmax = 30.0
    qmin = -1e-6
    indset = None
    print_opt = True
    cracks = True

    # Convex Optimisation

    fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
                                            printout=print_opt,
                                            find_inds=False,
                                            tol=0.01,
                                            objective='loadpath',
                                            indset=indset)

    form.to_json(file_initial)

    # Rollers or not
    if rollers:
        form = rollers_on_openings(form, xy_span = [[0.0,x_span],[0.0,y_span]], max_f = 5.0)

    # Sabouret Cracks
    form = create_cracks(form , dx =[[2.0, 8.0],[2.0, 8.0],[1.0, 1.0],[9.0, 9.0]], dy = [[9.0, 9.0],[1.0, 1.0],[2.0, 8.0],[2.0, 8.0]], type = ['top','top','top','top'], view = False)

    # Center Crack
    # form = create_cracks(form , dx =[[4.9, 5.1]], dy = [[4.9, 5.1]], type = ['top'], view = False)

    # Central Lines UB Crack
    # form = create_cracks(form , dx =[[5.0, 5.0],[0.0, 10.0]], dy = [[0.0, 10.0],[5.0, 5.0]], type = ['top','top'], view = False)

    # Irregular Crack Pattern
    # form = create_cracks(form , dx =[[2.0, 2.0],[2.5, 2.5],[7.0, 10.0],[5.0, 5.0]], dy = [[8.0, 10.0],[0.0, 2.5],[7.0, 7.0],[5.0,5.0]], type = ['bottom','bottom','bottom','top'], view = False)

    # plot_form(form).show()

    # form = FormDiagram.from_json(file_initial)
    # indset = form.attributes['indset']

    # Minimum Thrust

    solver = 'pyOpt-' + 'SLSQP'

    fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                        printout=print_opt,
                                        find_inds=True,
                                        indset=indset,
                                        tol=0.01,
                                        translation = translation,
                                        objective='lp_cracks',
                                        bmax = True,
                                        cracks = cracks,
                                        summary=print_opt)

    form.to_json(file_save)
    overview_forces(form)
    plot_form(form, show_q = False).show()
    reactions(form, plot = True)
