from compas_tna.diagrams import FormDiagram

from compas_thrust.utilities.symmetry import create_sym2
from compas_thrust.utilities.symmetry import replicate2
from compas_thrust.utilities.symmetry import create_sym
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import remove_feet
from compas_thrust.diagrams.form import oveview_forces
from compas_thrust.diagrams.form import adapt_objective

from compas_thrust.plotters.plotters import plot_form


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_comp.json'
    file_scaled = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_scaled.json'
    file_sym = '/Users/mricardo/compas_dev/compas_loadpath/data/freeform/SQ_sym.json'

    form = FormDiagram.from_json(file)

    # Scale to get an initial best-value for the objective with TNA

    form = adapt_objective(form, zrange = [3.0,8.0], objective = 'loadpath', plot = False, delete_face = True)
    form.to_json(file_scaled)
    oveview_forces(form)

    # Prepare Symmetrical part for optimisation with independents
    
    form = remove_feet(form, plot = False, openings = 12.0)

    # If form has 1 axis of symmetry

    # form = create_sym_2(form)
    # form.to_json(file_sym)

    # If form has 3 axis of Symmetry

    form = create_sym(form)
    form.to_json(file_sym)
    oveview_forces(form)