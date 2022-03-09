from compas_tno.algorithms import compute_reactions
from compas_tno.algorithms import xyz_from_q
from compas_tno.shapes import Shape
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.problems import save_geometry_at_iterations


def post_process_general(analysis):
    """Post processing of the optimisation.

    Parameters
    ----------
    analysis : Analysis
        The Analysis object

    Returns
    -------
    analysis : Analysis
        The Analysis object updated
    """

    form = analysis.form
    optimiser = analysis.optimiser
    shape = analysis.shape

    M = optimiser.M
    summary = optimiser.settings.get('summary', False)
    printout = optimiser.settings.get('printout', True)
    variables = optimiser.settings['variables']
    thickness_type = optimiser.settings.get('thickness_type', 'constant')
    features = optimiser.settings.get('features', [])
    save_iterations = optimiser.settings.get('save_iterations', False)

    fconstr = optimiser.fconstr
    # args = optimiser.args
    xopt = optimiser.xopt
    fopt = optimiser.fopt
    message = optimiser.message
    thk = M.thk

    i_uv = form.index_uv()
    # i_k = form.index_key()

    M.q = M.B.dot(xopt[:M.k])
    check = M.k

    if 'xyb' in variables:
        xyb = xopt[check:check + 2*M.nb]
        check = check + 2*M.nb
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in variables:
        zb = xopt[check: check + M.nb]
        check = check + M.nb
        M.X[M.fixed, [2]] = zb.flatten()
    if 't' in variables:
        thk = xopt[-1]
    if 'n' in variables:
        n = xopt[-1]
    if 'lambd' in variables:
        lambd = xopt[-1]
        M.P[:, [0]] = lambd * M.px0
        M.P[:, [1]] = lambd * M.py0
    if 'tub' in variables:
        tub = xopt[check:check + M.n]
        check = check + M.n
    if 'tlb' in variables:
        tlb = xopt[check:check + M.n]
        check = check + M.n
    if 'tub_reac' in M.variables:
        tub_reac = xopt[check: check + 2*M.nb]
        check = check + 2*M.nb
    # if 's' in variables:
    #     s = xopt[-1]

    g_final = fconstr(xopt, M)
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    i = 0
    for key in form.vertices():
        form.vertex_attribute(key, 'x', M.X[i, 0])
        form.vertex_attribute(key, 'y', M.X[i, 1])
        form.vertex_attribute(key, 'z', M.X[i, 2])
        form.vertex_attribute(key, 'px', M.P[i, 0])
        form.vertex_attribute(key, 'py', M.P[i, 1])
        form.vertex_attribute(key, 'pz', M.P[i, 2])
        i = i + 1

    for c, qi in enumerate(list(M.q.ravel())):
        u, v = i_uv[c]
        li = form.edge_length(u, v)
        form.edge_attribute((u, v), 'q', float(qi))
        form.edge_attribute((u, v), 'f', float(qi*li))

    form.attributes['loadpath'] = form.loadpath()
    compute_reactions(form)

    if 't' in variables:
        if shape.datashape['type'] == 'general':
            if thickness_type == 'constant':
                form.attributes['thk'] = thk
                shape.datashape['thk'] = thk
                shape.intrados = shape.middle.offset_mesh(n=thk/2, direction='down')
                shape.extrados = shape.middle.offset_mesh(n=thk/2, direction='up')
                apply_envelope_from_shape(form, shape)
            elif thickness_type == 'variable':
                t0 = shape.datashape['thk']
                thk = t0 * thk  # Consider that the thk for general shapes is a percentage of the thickness
                form.attributes['thk'] = thk
                shape.datashape['thk'] = thk
                if printout:
                    print('Optimum Value corresponds to a thickness of:', thk)
                shape.extrados, shape.intrados = shape.middle.offset_up_and_down(n=fopt)
                apply_envelope_from_shape(form, shape)
            elif thickness_type == 'intrados':
                form.attributes['thk'] = thk
                shape.datashape['thk'] = thk
                shape.middle = shape.intrados.offset_mesh(n=thk/2, direction='up')
                shape.extrados = shape.intrados.offset_mesh(n=thk, direction='up')
                form.envelope_from_shape(shape)
        else:
            form.attributes['thk'] = thk
            shape.datashape['thk'] = thk
            shape = Shape.from_library(shape.datashape)
            apply_envelope_from_shape(form, shape)  # Check if this is ok for adapted pattern
            i = 0
            for key in form.vertices():  # this resolve the problem due to the adapted pattern
                form.vertex_attribute(key, 'ub', float(M.ub[i]))
                form.vertex_attribute(key, 'lb', float(M.lb[i]))
                i += 1

    if 'update-envelope' in features:
        form.attributes['thk'] = thk
        shape.datashape['thk'] = thk
        shape = Shape.from_library(shape.datashape)
        apply_envelope_from_shape(form, shape)

    # if 's' in variables:
    #     s = -1 * fopt
    #     for key in form.vertices():
    #         ub = form.vertex_attribute(key, 'ub')
    #         lb = form.vertex_attribute(key, 'lb')
    #         form.vertex_attribute(key, 'ub', ub - s * (ub - lb))
    #         form.vertex_attribute(key, 'lb', lb + s * (ub - lb))

    if 'n' in variables:
        print('Value of N:', n)
        n = -1 * fopt
        shape.intrados = shape.intrados.offset_mesh(n=n, direction='up')
        shape.extrados = shape.extrados.offset_mesh(n=n, direction='down')
        apply_envelope_from_shape(form, shape)

    if 'tub' in variables:
        for i, key in enumerate(form.vertices()):
            zub = form.vertex_attribute(key, 'ub')
            form.vertex_attribute(key, 'tub', tub[i])
            form.vertex_attribute(key, 'ub', zub + tub[i])

    if 'tlb' in variables:
        for i, key in enumerate(form.vertices()):
            zub = form.vertex_attribute(key, 'lb')
            form.vertex_attribute(key, 'tlb', tlb[i])
            form.vertex_attribute(key, 'lb', zub - tlb[i])

    if 'tub_reac' in variables:
        for i, key in enumerate(form.vertices_where({'is_fixed': True})):
            print(i, key, tub_reac)
            form.vertex_attribute(key, 'tub_reac', [tub_reac[i], tub_reac[i + M.nb]])

    analysis.form = form
    analysis.optimiser = optimiser
    analysis.shape = shape

    if save_iterations:
        save_geometry_at_iterations(form, optimiser, force=True)

    if printout or summary:
        print('\n' + '-' * 50)
        print('Solution  :', message)
        try:
            print('q range : {0:.3f} : {1:.3f}'.format(min(M.q), max(M.q)))
        except BaseException:
            print('q range : {0:.3f} : {1:.3f}'.format(min(M.q.flatten()), max(M.q.flatten())))
        print('zb range  : {0:.3f} : {1:.3f}'.format(min(M.X[M.fixed, [2]]), max(M.X[M.fixed, [2]])))
        print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
        print('fopt      : {0:.3f}'.format(fopt))
        print('-' * 50 + '\n')

    return analysis
