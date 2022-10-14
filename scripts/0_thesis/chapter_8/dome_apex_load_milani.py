
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.viewers import Viewer
from numpy import zeros

r = 1.21
thk = 0.12

form = FormDiagram.create_circular_radial_form(discretisation=[16, 20], radius=r)
shape = Shape.create_dome(discretisation=[16, 20], radius=r, thk=thk)
shape.ro = 18.0

key = 0

load_direction = zeros([form.number_of_vertices(), 1])
load_direction[key] = - 1.0

analysis = Analysis.create_max_load_analysis(form, shape,
                                             load_direction=load_direction,
                                             max_lambd=1000.0)
analysis.optimiser.set_constraints(['funicular', 'envelope', 'reac_bounds'])
analysis.apply_selfweight()
analysis.apply_envelope()
lumped_swt = form.lumped_swt()
print('Weight Shape:', shape.compute_selfweight())
print('SWT', lumped_swt)
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# lumped_swt = form.lumped_swt()
pload = analysis.optimiser.fopt

print('Opt. Val:', pload)
print('Weight Shape:', shape.compute_selfweight())
print('SWT', lumped_swt)
print('Pmax/W:', pload/lumped_swt)

view = Viewer(form, shape)
view.scale_edge_thickness(10.0)
view.draw_shape()
view.draw_thrust()
view.draw_cracks()
view.draw_reactions(extend_reactions=True)
view.show()
