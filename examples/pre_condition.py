from compas_tna.diagrams import FormDiagram

from compas_thrust.utilities.symmetry import create_sym2
from compas_thrust.utilities.symmetry import replicate2
from compas_thrust.utilities.symmetry import create_sym
from compas_thrust.utilities.symmetry import replicate

from compas.utilities import geometric_key

from compas_thrust.diagrams.form import remove_feet
from compas_thrust.diagrams.form import oveview_forces
from compas_thrust.diagrams.form import adapt_objective
from compas_thrust.diagrams.form import adapt_tna

from compas_thrust.plotters.plotters import plot_form

from compas_viewers.meshviewer import MeshViewer


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file_complete = '/Users/mricardo/compas_dev/me/minmax/radial/mixed_05_complete.json'
    # file = '/Users/mricardo/compas_dev/me/loadpath/prototype/prot_complete.json'
    file_scaled = '/Users/mricardo/compas_dev/me/minmax/radial/mixed_05_scaled.json'
    # file_scaled = '/Users/mricardo/compas_dev/me/loadpath/prototype/prot_scaled.json'
    file_sym = '/Users/mricardo/compas_dev/me/minmax/radial/mixed_05_init.json'
    # file_real = '/Users/mricardo/compas_dev/me/minmax/radial/mixed_05_sym.json'

    form = FormDiagram.from_json(file_complete)
    plot_form(form).show()

    # Scale to get an initial best-value for the objective with TNA

    # form = adapt_objective(form, zrange = [3.0,9.0], kmax = 5000, objective = 'target', plot = True, delete_face = True)
    # form.to_json(file_scaled)
    # plot_form(form).show()

    # Prepare Symmetrical part for optimisation with independents
    
    # form = remove_feet(form, plot = True)
    # oveview_forces(form)
    # form.to_json(file_scaled)

    # for key in form.vertices():
    #     print(form.get_vertex_attribute(key, 'target'))

    # If form has 1 axis of symmetry

    # form = create_sym_2(form)
    # form.to_json(file_sym)

    # If form has 3 axis of Symmetry

    form = create_sym(form)
    form.to_json(file_sym)
    plot_form(form).show()
    oveview_forces(form)

    # Viewer

    viewer = MeshViewer()
    viewer.mesh = form
    viewer.show()

    # Update Qs

    # qs = {}

    # for u,v in form.edges():
    #     gkey = geometric_key(form.edge_midpoint(u,v)[:2] + [0])
    #     qs[gkey] = form.get_edge_attribute((u,v), 'q')
    
    # print(qs)

    # form_ = FormDiagram.from_json(file_real)

    # plot_form(form_).show()

    # for key, attr in form_.vertices(True):
    #     attr['z'] = 0.0

    # for u,v in form_.edges():
    #     gkey = geometric_key(form_.edge_midpoint(u,v)[:2] + [0])
    #     form_.set_edge_attribute((u,v), name='q', value=qs[gkey])

    # plot_form(form_).show()
    # form_.to_json(file_real)

    #Replicate

    form_ = replicate(form, file_complete, plot = True)
    # form_ = replicate2(form, file)
    form_.to_json(file_complete)

    # Viewer

    viewer = MeshViewer()
    viewer.mesh = form_
    viewer.show()