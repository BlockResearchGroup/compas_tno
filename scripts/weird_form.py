from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import initialise_form
from compas_tno.algorithms import initialise_problem
from compas_tno.viewers import view_thrust
from compas_tno.algorithms import zlq_from_qid

address = '/Users/mricardo/compas_dev/me/reformulation/weird_form.json'
form = FormDiagram.from_json(address)
args = initialise_problem(form)
# form = initialise_form(form)
plot_form(form).show()
# view_thrust(form).show()

q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

qid = q[ind]
z_, _, _, _ = zlq_from_qid(qid, args)

differences = z - z_
print(max(differences))
print(min(differences))
