from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form

# ------------------------------------------------
# --------- CREATING ARCH FORM DIAGRAM -----------
# ------------------------------------------------

data = {
    'type': 'arch',
    'H': 1.0,
    'L': 2.0,
    'total_nodes': 11,
    'x0': 0.0,
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()

# ------------------------------------------------
# --------- CREATING ORTHO FORM DIAGRAM ----------
# ------------------------------------------------

data = {
    'type': 'ortho',
    'xy_span': [[0,10],[0,10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()

# ------------------------------------------------
# --------- CREATING CROSS FORM DIAGRAM ----------
# ------------------------------------------------

data = {
    'type': 'cross_fd',
    'xy_span': [[0,10],[0,10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()

# ------------------------------------------------
# --------- CREATING FAN FORM DIAGRAM ------------
# ------------------------------------------------

data = {
    'type': 'fan_fd',
    'xy_span': [[0,10],[0,10]],
    'discretisation': [10, 10],
    'fix': 'corners',
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()

# ------------------------------------------------
# --------- CREATING RADIAL FORM DIAGRAM ---------
# ------------------------------------------------

data = {
    'type': 'radial_fd',
    'D': 3.0,
    'center': [5.0, 5.0],
    'radius': 5.0,
    'discretisation': [8, 20],
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()
# ------------------------------------------------
# --------- CREATING RADIAL SPACED DIAGRAM ---------
# ------------------------------------------------

data = {
    'type': 'radial_spaced_fd',
    'D': 3.0,
    'center': [5.0, 5.0],
    'radius': 5.0,
    'discretisation': [8, 20],
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()


# ------------------------------------------------
# --------- CREATING SPIRAL FORM DIAGRAM ---------
# ------------------------------------------------

data = {
    'type': 'spiral_fd',
    'D': 3.0,
    'center': [5.0, 5.0],
    'radius': 5.0,
    'discretisation': [8, 20],
    'r_oculus': 0.0,
}

form = FormDiagram.from_library(data)
print(form)
form.overview_forces()
plot_form(form, show_q=False).show()
