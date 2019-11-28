from compas_tna.diagrams import FormDiagram

from compas_thrust.utilities.constraints import circular_heights
from compas_thrust.utilities.constraints import circular_joints
from compas_thrust.utilities.constraints import set_pointed_vault_heights
from compas_thrust.utilities.constraints import set_pointed_vault_heights_ub
from compas_thrust.utilities.constraints import set_pointed_vault_heights_lb
from compas_thrust.utilities.constraints import set_cross_vault_heights
from compas_thrust.utilities.constraints import set_cross_vault_heights_ub
from compas_thrust.utilities.constraints import set_cross_vault_heights_lb
from compas_thrust.utilities.constraints import create_cracks
from compas_thrust.utilities.constraints import rollers_on_openings


from compas_thrust.diagrams.form import overview_forces
from compas_thrust.diagrams.form import create_arch
from compas_thrust.diagrams.form import create_cross_form
from compas_thrust.diagrams.form import create_fan_form
from compas_thrust.diagrams.form import _form
from compas_thrust.algorithms import optimise_general
from compas_thrust.algorithms import optimise_convex

from compas_thrust.plotters.plotters import plot_form_xz
from compas_thrust.plotters.plotters import plot_form_joints
from compas_thrust.plotters.plotters import plot_form
from numpy import array

from compas.geometry import intersection_segment_segment_xy
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import intersection_line_line_xy

from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import array
from numpy.linalg import matrix_rank
from numpy import transpose
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import scale_vector_xy
from compas.geometry import distance_point_point_xy
from compas_thrust.algorithms import zlq_from_qid
from compas_thrust.algorithms import q_from_qid
from scipy.sparse import diags

from compas_viewers.multimeshviewer import MultiMeshViewer

from compas_thrust.algorithms.problems import initialise_problem

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file_save = '/Users/mricardo/compas_dev/me/minmax/ortho/square/D_20_min_thrust+cracks.json'
    file_initial = '/Users/mricardo/compas_dev/me/minmax/ortho/square/D_20_min_loadpath.json'


    # Create Vault from one of the patterns Fan/Grid and set constraints to be cross-vault
    
    x_span = 10.0
    y_span = 10.0
    divisions = 4
    form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
    form = set_cross_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], thk = 0.5, b = 5.0, set_heights=False, ub_lb = True, update_loads = True)
    plot_form(form).show()
    # form = create_cracks(form , dx =[[4.9, 5.1]], dy = [[4.9, 5.1]], type = ['top'], view = False)
    # form = create_cracks(form , dx =[[2.0, 2.0],[2.0, 2.0],[8.0, 8.0],[8.0, 8.0],[5.0, 5.0]], dy = [[2.0, 2.0],[8.0, 8.0],[8.0, 8.0],[2.0, 2.0],[5.0, 5.0]], type = ['bottom','bottom','bottom','bottom','top'], view = False)
    # form = create_cracks(form , dx =[[5.0, 5.0],[0.0, 10.0]], dy = [[0.0, 10.0],[5.0, 5.0]], type = ['top','top'], view = False)
    # form = create_cracks(form , dx =[[2.0, 2.0],[2.5, 2.5],[7.0, 10.0],[5.0, 5.0]], dy = [[8.0, 10.0],[0.0, 2.5],[7.0, 7.0],[5.0,5.0]], type = ['bottom','bottom','bottom','top'], view = False)
    form = rollers_on_openings(form)

    # Initial parameters

    translation = form.attributes['tmax']
    bounds_width = 5.0
    use_bounds = False
    qmax = 100
    qmin = -1e-6
    indset = None
    print_opt = True
    cracks = True

    # Convex Optimisation

    # fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
    #                                         printout=print_opt,
    #                                         find_inds=True,
    #                                         tol=0.01,
    #                                         objective='loadpath',
    #                                         indset=indset)

    # form.to_json(file_initial)

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
                                        use_bounds = use_bounds,
                                        bounds_width = bounds_width,
                                        objective='max',
                                        bmax = True,
                                        cracks = cracks,
                                        summary=print_opt)

    form.to_json(file_save)

    # Visualisation
    
    # W.I.P.