from compas_tno.algorithms import reactions
from compas_tno.algorithms import zq_from_qid

from compas.utilities import geometric_key

__all__ = [
    'post_process_analysis'
]


def post_process_analysis(analysis):

    form = analysis.form
    optimiser = analysis.optimiser

    plot = optimiser.data.get('plot', False)
    summary = optimiser.data.get('summary', False)
    printout = optimiser.data.get('printout', True)
    variables = optimiser.data['variables']

    fconstr = optimiser.fconstr
    args = optimiser.args
    xopt = optimiser.xopt
    fopt = optimiser.fopt
    message = optimiser.message

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k = args[:23]

    i_uv = form.index_uv()
    i_k = form.index_key()

    if 'ind' in variables:
        q[ind] = xopt[:k].reshape(-1, 1)
    else:
        q = xopt[:len(q)].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)
    if 't' in variables:
        thk = xopt[-1]
    if 's' in variables:
        s = xopt[-1]

    g_final = fconstr(xopt, *args)
    z, q = zq_from_qid(q[ind], args)

    gkeys = []
    for i in ind:
        u, v = i_uv[i]
        gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
    form.attributes['indset'] = gkeys

    for i in range(form.number_of_vertices()):
        key = i_k[i]
        form.vertex_attribute(key, 'z', float(z[i]))

    for c, qi in enumerate(list(q.ravel())):
        u, v = i_uv[c]
        form.edge_attribute((u, v), 'q', float(qi))

    form.attributes['loadpath'] = form.loadpath()
    reactions(form, plot=plot)

    if 't' in variables:
        form.attributes['thk'] = thk

    if 's' in variables:
        s = -1 * fopt
        for key in form.vertices():
            ub = form.vertex_attribute(key, 'ub')
            lb = form.vertex_attribute(key, 'lb')
            form.vertex_attribute(key, 'ub', ub - s * (ub - lb))
            form.vertex_attribute(key, 'lb', lb + s * (ub - lb))

    analysis.form = form
    analysis.optimiser = optimiser

    if printout or summary:
        print('\n' + '-' * 50)
        print('Solution  :', message)
        print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
        print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
        print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
        print('fopt      : {0:.3f}'.format(fopt))
        print('-' * 50 + '\n')

    return analysis
