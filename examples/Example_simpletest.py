from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_form_xz
from compas_tno.shapes import Shape

file_address = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form = FormDiagram.from_json(file_address)

# Basic parameters

H = 0.5
L = 2.0
thk = 0.2
radius = 5.0
discretisation = 20
b = 0.5
t = 5.0
type_structure = 'arch'
type_formdiagram = 'arch'

# ----------------------- 1. Create Arch shape ---------------------------

data_shape = {
    'type': type_structure,
    'H': H,
    'L': L,
    'thk': thk,
    'discretisation': discretisation,
    'b': b,
    't': t,
    'x0': 0.0
}

arch = Shape.from_library(data_shape)
swt = arch.compute_selfweight()
print('Arch created!')
# view_shapes(arch).show()

for key in form.vertices():
    pz = form.vertex_attribute(key, 'pz')
    print(pz)

for key in form.vertices_where({'is_fixed': True}):
    print(key)
    pz = form.vertex_attribute(key, 'pz')
    print(pz)

plot_form_xz(form, arch, show_q=False, plot_reactions=True, fix_width=True, max_width=5, radius=0.02).show()
