import json
import compas_tno
from numpy import array

from compas_tno.algorithms import xyz_from_xopt
from compas_tno.algorithms import reciprocal_from_form


def callback_save_json(xopt, *args, **kwargs):
    """Save the variables in the ``output.json`` file created.

    Parameters
    ----------
    xopt : array
        The variables in one iteration of the optimisation
    """

    DATA_FILENAME = compas_tno.get('output.json')

    with open(DATA_FILENAME, mode='r', encoding='utf-8') as f:
        data = json.load(f)

    i = len(data['iterations'])
    data['iterations'][i] = list(xopt.flatten())

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return


def callback_create_json():
    """Create a ``output.json`` to store the iterations of an optimisation"""

    data = {'iterations': {}}

    DATA_FILENAME = compas_tno.get('output.json')

    with open(DATA_FILENAME, mode='w', encoding='utf-8') as f:
        json.dump(data, f)

    return


def save_geometry_at_iterations(form, optimiser, force=False):
    """Save the geometry of the form (and force) during iterations of the optimisation. Works only with SLSQP and IPOPT solvers.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The Form Diagram of the problem
    optimiser : :class:`~compas_tno.optimisers.Optimiser`
        The Optimiser with the numerical information about the problem
    force : bool, optional
        Whether or not the force diagram coordinates should also be saved, by default False

    Notes
    -----
        Xforce
            Address where the geometry of the structure at each step is saved ``Xform.json``
        Xforce
            Address where the geometry of the force diagram at each step is saves ``Xforce.json``.

    """

    M = optimiser.M  # matrices of the problem

    file_qs = compas_tno.get('output.json')
    file_Xform = compas_tno.get('Xform.json')
    file_Xforce = None

    if force:
        file_Xforce = compas_tno.get('Xforce.json')
        Xforce = {}

    with open(file_qs, mode='r', encoding='utf-8') as f:
        data = json.load(f)

    Xform = {}

    iterations = len(data['iterations'])

    for i in range(iterations):
        xopt_i = array(data['iterations'][str(i)]).reshape(-1, 1)
        M = xyz_from_xopt(xopt_i, M)
        Xform_i = M.X.tolist()
        Xform[str(i)] = Xform_i

        if force:
            j = 0
            for key in form.vertices():
                form.vertex_attributes(key, 'xyz', M.X[j].tolist())
                j += 1
            k = 0
            for edge in form.edges_where({'_is_edge': True}):
                form.edge_attribute(edge, 'q', M.q.flatten()[k])
                k += 1

            force = reciprocal_from_form(form)
            Xforce_i = force.vertices_attributes('xyz')
            Xforce[str(i)] = Xforce_i

    with open(file_Xform, mode='w', encoding='utf-8') as f:
        json.dump(Xform, f)

    print('Form Geometry saved @:', file_Xform)

    if force:
        with open(file_Xforce, mode='w', encoding='utf-8') as f:
            json.dump(Xforce, f)

        print('Force Geometry saved @:', file_Xforce)

    return file_Xform, file_Xforce
