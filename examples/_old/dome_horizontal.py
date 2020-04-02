from compas_tna.diagrams import FormDiagram

from compas_tno.diagrams.form import overview_forces
from compas_tno.diagrams.form import create_cross_form
from compas_tno.diagrams.form import create_fan_form
from compas_tno.diagrams.form import create_dome_form
from compas_tno.diagrams.form import delete_boundary_edges

from compas_tno.utilities.constraints import set_pavillion_vault_heights
from compas_tno.utilities.constraints import set_dome_heights

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import optimise_general
from compas_tno.algorithms import optimise_convex

from compas.datastructures import mesh_quads_to_triangles

from compas_tno.plotters import plot_form

import math


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    # Try with 'fan_fd' and 'cross_fd' and for the objective change 'min' and 'max'
    type_fd = 'radial'
    objective = 'feasibility'
    thck = 0.30

    # Create Vault from one of the patterns Fan/Grid with the dimensions

    xc = 5.0
    yc = 5.0
    radius = 5.0
    n_radial = 8
    n_spikes = 16
    rx = 0.01

    form = create_dome_form(center = [xc, yc], radius = radius, n_radial = n_radial, n_spikes = n_spikes, r_oculus= 0.0)
    form = set_dome_heights(form, center = [xc, yc], radius = radius, thck = thck)
    form = delete_boundary_edges(form)

    PATH = '/Users/mricardo/compas_dev/me/minmax/dome/hor-loads/r=' + str(int(radius)) + '/' + type_fd + '_discr_'+ str(n_radial) + '_' + str(n_spikes)
    file_initial = PATH + '_lp.json'
    file_save = PATH + '_' + objective + '_t=' + str(int(thck*100)) + '_px=' + str(rx) + '.json'

    translation = True
    qmax = 200
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
    form = FormDiagram.from_json(file_initial)
    plot_form(form, show_q = False, fix_width=True).show()
    # indset = form.attributes['indset']

    # Add Horizontal Loads

    mesh_quads_to_triangles(form)
    plot_form(form, show_q = False, fix_width=True).show()

    # for u,v in form.edges():
    #     qi = form.edge_attribute((u,v),'q')
    #     if qi == 1.0 or qi == None:
    #         form.edge_attribute((u,v),'q',value=0.0)

    pxt = 0
    for key in form.vertices():
        form.vertex[key]['px'] = form.vertex[key]['pz'] * rx
        pxt += form.vertex[key]['px']
    print('Load x-direction - pxt = {0}'.format(pxt))

    # Maximum or Minimum Thrusts

    solver = 'pyOpt-' + 'SLSQP'
    # solver = 'slsqp'

    fopt, qopt, zbopt, exitflag = optimise_general(form,  qmax=qmax, solver=solver,
                                        printout=print_opt,
                                        find_inds=True,
                                        indset=indset,
                                        translation = translation,
                                        objective=objective,
                                        bmax = True,
                                        summary=print_opt)

    print('Form Saved to: ',file_save)
    form.to_json(file_save)
    overview_forces(form)
    plot_form(form, show_q = False).show()
    reactions(form)

    for key in form.vertices_where({'is_fixed': True}):
        rx = round(form.vertex_attribute(key, '_rx'),3)
        ry = round(form.vertex_attribute(key, '_ry'),3)
        zb = round(form.vertex_attribute(key,'z'),3)
        print('Reaction on Corner {0}: rx: {1:.3f}/ ry: {2:.3f}/ r: {3:.3f}'.format(key, rx, ry, math.sqrt(rx**2 + ry**2)))
        break
    qmax = round(max(qopt).item(),3)
    fopt = round(fopt,3)
    print('zb: {0:.3f}'.format(zb))
    print('fopt: {0:.3f}'.format(fopt))
    print('qmax: {0:.3f}'.format(max(qopt).item()))
