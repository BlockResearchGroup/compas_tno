from compas_tna.diagrams import ForceDiagram

title = 'pointed_crossvault_cross_fd_t=50_min_force'
force_address = '/Users/mricardo/compas_dev/compas_tno/data/rqe/' + title + '.json'
force = ForceDiagram.from_json(force_address)

force.draw(layer=title)