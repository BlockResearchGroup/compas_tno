from compas_tna.diagrams import FormDiagram

from compas_thrust.utilities.constraints import circular_heights
from compas_thrust.utilities.constraints import circular_joints
from compas_thrust.utilities.constraints import set_pointed_vault_heights
from compas_thrust.utilities.constraints import set_pointed_vault_heights_ub
from compas_thrust.utilities.constraints import set_pointed_vault_heights_lb
from compas_thrust.utilities.constraints import set_cross_vault_heights
from compas_thrust.utilities.constraints import set_cross_vault_heights_ub
from compas_thrust.utilities.constraints import set_cross_vault_heights_lb

from compas_thrust.diagrams.form import overview_forces
from compas_thrust.diagrams.form import create_arch
from compas_thrust.diagrams.form import create_cross_form
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

    # file_pattern = '/Users/mricardo/compas_dev/me/loadpath/freeform/test.json'
    # '/Users/mricardo/compas_dev/me/loadpath/corner/discretize/01_08_complete_paper.json'
    file_save = '/Users/mricardo/compas_dev/me/loadpath/corner/pointed/rounded.json'
    # form = FormDiagram.from_json(file_pattern)
    # for key in form.edges():
    #     form.set_edge_attribute(key, 'q', value = 1.0)
    # for key in form.vertices():
    #     form.set_vertex_attribute(key, 'z', value = 0.0)
    # form.to_json('/Users/mricardo/Documents/MATLAB/optimisation/discretize/form1.json')
    # plot_form(form, show_q=False, simple=True).show()
    x_span = 10.0
    y_span = 5.0
    divisions = 50
    form = create_cross_form(xy_span = [[0.0,x_span],[0.0,y_span]], division=divisions)
    # form = set_cross_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], set_heights=True)
    form = set_cross_vault_heights_lb(form, xy_span = [[0.0,x_span],[0.0,y_span]], set_heights=True, thk = 0.4)

    # plot_form(form, show_q=False, simple=True, heights=True).show()
    # form.to_json('/Users/mricardo/Documents/MATLAB/optimisation/discretize/form2.json')
    # print('Number of Edges:',form.number_of_edges())
    # form = _form(form)
    # form.to_json('/Users/mricardo/Documents/MATLAB/optimisation/discretize/form2.json')
    # plot_form(form, show_q=False, simple=True).show()
    # form = _form(form)
    # plot_form(form, show_q=False, simple=True).show()
    # print(form)
    # plot_form(form, show_q=True).show()
    # form = set_pointed_vault_heights_lb(form, xy_span = [[0.0,x_span],[0.0,y_span]], hc=8.0, he=[6.0,6.0,6.0,6.0], set_heights=True, thk = 0.6)
    
    # form = set_pointed_vault_heights_lb(form, xy_span = [[0.0,x_span],[0.0,y_span]], hc=6.0, he=[7.0,7.0,7.0,7.0], hm=[8.0,8.0,8.0,8.0], set_heights=True, thk = 0.6)
    # form = set_pointed_vault_heights(form, xy_span = [[0.0,x_span],[0.0,y_span]], hc=7.0, set_heights=True, thk = 0.6)

    # Initial parameters

    translation = None
    bounds_width = 5.0
    use_bounds = False
    qmax = 20
    indset = None
    print_opt = True

    # Optimisation

    # fopt, qopt, zbopt, exitflag = optimise_convex(form,  qmax=qmax,
    #                                         printout=print_opt,
    #                                         find_inds=True,
    #                                         tol=0.01,
    #                                         objective='loadpath',
    #                                         indset=indset)

    form.to_json(file_save)

    # Visualisation

    # from compas.datastructures import Mesh

    # from compas_viewers import VtkViewer
    # import compas

    # # datastructure = Mesh.from_obj(compas.get('quadmesh.obj'))

    # viewer = VtkViewer(datastructure=form)
    # viewer.setup()
    # viewer.start()



    