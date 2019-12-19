from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tno.plotters.plotters import plot_form
from compas_tno.plotters.plotters import plot_force
from compas_tno.algorithms import initialize_problem
from compas_tno.diagrams.form import overview_forces
from compas_tno.algorithms import update_tna



# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    i_s = []
    n_s = []
    lp_s = []
    i = 5
    j = 2

    file = '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/0'+ str(j) +'_0'+ str(i) +'_complete_paper.json'
    form = FormDiagram.from_json(file)
    args = initialize_problem(form)
    plot_form(form, show_q=False, max_width=3.0).show()
    force = ForceDiagram.from_formdiagram(form)
    form, force = update_tna(form, delete_face=True)
    plot_force(force, form, color_inds=True).show()

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args

    print(ind)
    # q[ind] = [ change as you want]
    # q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
    # Assign the vector q to you mesh as for.set_edge_attribute(key, 'q', value = q[i])
