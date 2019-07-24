from compas_tna.diagrams import FormDiagram
from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.plotters.plotters import plot_form

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    file = '/Users/mricardo/compas_dev/me/minmax/barrel/2D_min.json'
    form = FormDiagram.from_json(file)
    plot_form(form)
    reactions(form, plot=True)
    form.to_json(file)