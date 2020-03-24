from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights

from compas_tno.algorithms.equilibrium import reactions

from compas.datastructures import mesh_quads_to_triangles

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas_viewers.meshviewer import MeshViewer

from compas_tno.plotters import plot_form

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'cross_fd'
    objective = 'constr_lp'
    objective_load = 'loadpath'
    thck = 0.50

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    x_span = 10.0
    y_span = 10.0
    px = 0.1

    if type_fd == 'cross_fd':
        divisions = 12
        xy_span = [[0.0, x_span], [0.0, y_span]]
        form = create_cross_form(xy_span=[[0.0, x_span], [0.0, y_span]], division=divisions, fix='all', px = px)  # FIX ALL NODES ON BOUNDARIES
    if type_fd == 'fan_fd':
        divisions = 16
        xy_span = [[0.0, x_span], [0.0, y_span]]
        form = create_fan_form(xy_span=[[0.0, x_span], [0.0, y_span]], division=divisions, fix='all')  # FIX ALL NODES ON BOUNDARIES

    form = delete_boundary_edges(form)
    mesh_quads_to_triangles(form)

    PATH = '/Users/mricardo/compas_dev/me/minmax/pavillion/hor-loads/' + type_fd + '/' + type_fd + '_discr_' + str(divisions)

    # file_initial = PATH + '_lp.json'
    file_initial = PATH + '_' + objective_load + '_t=' + str(int(thck*100)) + '_px=' + str(px) + '.json'

    # thk = 1.0
    file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '_px=' + str(px) + '.json'

    form = FormDiagram.from_json(file_initial)

    # One run for loadpath

    # Initial parameters

    translation = False
    qmax = 400
    qmin = -1e-6
    indset = None
    print_opt = True

    # Convex Optimisation to find good starting point. Save the starting point, and can load it later if wanted

    form = set_pavillion_vault_heights(form, xy_span=[[0.0, x_span], [0.0, y_span]], thk=thck, b=5.0, t=0.0, set_heights=False, ub_lb=True, update_loads=True)
    plot_form(form, show_q=False).show()

    # Maximum or Minimum Thrusts

    solver = 'pyOpt-' + 'SLSQP'
    # solver = 'slsqp'

    fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                                   printout=print_opt,
                                                   find_inds=True,
                                                   indset=indset,
                                                   translation=translation,
                                                   objective=objective,
                                                   bmax=True,
                                                   summary=print_opt)

    print('File Saved:', file_save)
    form.to_json(file_save)
    overview_forces(form)
    # plot_form(form, show_q = False).show()
    reactions(form)

    for key in form.vertices_where({'is_fixed': True}):
        rx = round(form.vertex_attribute(key, 'rx'), 3)
        ry = round(form.vertex_attribute(key, 'ry'), 3)
        zb = round(form.vertex_attribute(key, 'z'), 3)
        print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
        break
    qmax = round(max(qopt).item(), 3)
    fopt = round(fopt, 3)
    print('zb: {0:.3f}'.format(zb))
    print('fopt: {0:.3f}'.format(fopt))
    print('qmax: {0:.3f}'.format(max(qopt).item()))

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()