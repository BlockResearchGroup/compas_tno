from compas_tno.algorithms import reactions
from compas_tno.algorithms import zq_from_qid
from compas_tno.shapes import Shape

from compas.utilities import geometric_key

__all__ = [
    'post_process_analysis'
]


def post_process_analysis(analysis):

    form = analysis.form
    optimiser = analysis.optimiser
    shape = analysis.shape

    plot = optimiser.data.get('plot', False)
    summary = optimiser.data.get('summary', False)
    printout = optimiser.data.get('printout', True)
    variables = optimiser.data['variables']
    thickness_type = optimiser.data.get('thickness_type', 'constant')

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
        if shape.data['type'] == 'general':
            if thickness_type == 'constant':
                form.attributes['thk'] = thk
                shape.data['thk'] = thk
                shape.intrados = shape.middle.offset_mesh(n=thk/2, direction='down')
                shape.extrados = shape.middle.offset_mesh(n=thk/2, direction='up')
                form.envelope_from_shape(shape)
            elif thickness_type == 'variable':
                t0 = shape.data['thk']
                thk = t0 * thk  # Consider that the thk for general shapes is a percentage of the thickness
                form.attributes['thk'] = thk
                shape.data['thk'] = thk
                if printout:
                    print('Optimum Value corresponds to a thickness of:', thk)
                shape.extrados, shape.intrados = shape.middle.offset_up_and_down(n=fopt)
                form.envelope_from_shape(shape)
            elif thickness_type == 'intrados':
                form.attributes['thk'] = thk
                shape.data['thk'] = thk
                shape.middle = shape.intrados.offset_mesh(n=thk/2, direction='up')
                shape.extrados = shape.intrados.offset_mesh(n=thk, direction='up')
                form.envelope_from_shape(shape)
        else:
            form.attributes['thk'] = thk
            shape.data['thk'] = thk
            shape = Shape.from_library(shape.data)
            form.envelope_from_shape(shape)

    if 's' in variables:
        s = -1 * fopt
        for key in form.vertices():
            ub = form.vertex_attribute(key, 'ub')
            lb = form.vertex_attribute(key, 'lb')
            form.vertex_attribute(key, 'ub', ub - s * (ub - lb))
            form.vertex_attribute(key, 'lb', lb + s * (ub - lb))

    if 'n' in variables:
        n = -1 * fopt
        shape.intrados = shape.intrados.offset_mesh(n=n, direction='up')
        shape.extrados = shape.extrados.offset_mesh(n=n, direction='down')
        form.envelope_from_shape(shape)

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
