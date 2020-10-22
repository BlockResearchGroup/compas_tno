from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.viewers.thrust import view_thrust
from compas_tno.plotters import plot_form

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR DOME --------------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 16]

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'center': [5.0, 5.0],
    'radius': radius,
    't' : 1.0
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()

print('Vault geometry created!')

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=False).show()


# --------------------- 3. Load the Form Diagram according to shape ---------------------

form.selfweight_from_shape(dome)
form.initialise_loadpath()

plot_form(form, show_q=False).show()
view_thrust(form).show()
