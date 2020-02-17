from compas_tno.diagrams.form import create_dome_flower
from compas_tno.diagrams.form import create_dome_form
from compas_tno.diagrams.form import create_dome_form_spaced
from compas_tno.plotters.plotters import plot_form
from compas_tno.algorithms.problems import initialise_problem
from compas_tno.diagrams.form import delete_boundary_edges
from compas_tna.diagrams import FormDiagram
from compas_tno.plotters.plotters import plot_form_xz
from compas_plotters import MeshPlotter

# xc = 5.0
# yc = 5.0
# radius = 5.0
# n_radial = 8
# n_spikes = 20

# form = create_dome_flower(center=[xc, yc], radius=radius, n_radial=n_radial, n_spikes=n_spikes, r_oculus=0.0)
# form = delete_boundary_edges(form)
# args = initialise_problem(form)
# plot_form(form, show_q = False, fix_width= True).show()


type_vault = 'dome'
type_fd = 'radial_spaced'
objective = 'max'
thck = 0.50
R = 5.0

xc = 5.0
yc = 5.0
radius = 5.0
n_radial = 8
n_spikes = 20
tol = 0.01

PATH = '/Users/mricardo/compas_dev/me/minmax/' + type_vault + '/' + type_fd + '/' + type_fd + '_discr_' + str(n_radial) + '_' + str(n_spikes)
try:
    file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*1000))) + '.json'
    print(file_open)
    form = FormDiagram.from_json(file_open)
except:
    file_open = PATH + '_' + objective + '_t=' + str(int(round(thck*100))) + '.json'
    print(file_open)
    form = FormDiagram.from_json(file_open)
print(form.number_of_edges())
# plot_form(form, show_q = False, fix_width= False).show()
form.attributes['Re'] = R + thck/2
form.attributes['Ri'] = R - thck/2
ind_total = 0
for key in form.edges_where({'is_ind': True}):
    ind_total += 1
print('Total of Independents: ', ind_total)
# plot_form(form, fix_width=False, show_edgeuv=False, heights= False, show_q=False).show()
# plotter = MeshPlotter(form, figsize=(10, 10))
# plotter.show()
# form.plot()
plot_form_xz(form, radius = 0.06, simple=True, fix_width = True, max_width=1.5, heights=True, show_q=False, thk = thck, plot_reactions=True, yrange = [R-tol,R+tol]).show()
