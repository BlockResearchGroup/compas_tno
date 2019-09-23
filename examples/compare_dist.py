
from compas_tna.diagrams import FormDiagram

from compas_thrust.algorithms.ind_based import optimise_single

from compas_thrust.algorithms.equilibrium import reactions
from compas_thrust.algorithms.equilibrium import horizontal_check

from compas_thrust.utilities.constraints import check_constraints
from compas_thrust.diagrams.form import overview_forces
from compas_thrust.utilities.symmetry import replicate

from compas_thrust.diagrams.form import _form

from compas.geometry import closest_point_in_cloud
from compas.geometry import distance_point_point
from compas.utilities import geometric_key
from compas.utilities import i_to_white
from compas_plotters import MeshPlotter

from compas_thrust.plotters.plotters import plot_form
import compas_pattern

from copy import deepcopy
from numpy import array
from numpy import argmin


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    averages = []
    maxs = []
    i_s = []
    
    for i in range(2,9):
        j = 2
        # file = '/Users/mricardo/compas_dev/me/discretize/0'+str(j)+'_0'+str(i)+'_complete.json'
        file_compare = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_09_complete.json'
        # file_compare = '/Users/mricardo/compas_dev/me/discretize/02_09_complete.json'
        file_complete = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_complete.json'
        file_dist = '/Users/mricardo/compas_dev/me/loadpath/Fix/discretize/0'+str(j)+'_0'+str(i)+'_dist.json'
        
        form = FormDiagram.from_json(file_dist)
        # form = FormDiagram.from_json(file)
        # formbase = FormDiagram.from_json(file_compare)
        # formbase_0 = deepcopy(formbase)

        # form.update_default_vertex_attributes({'dist': 0.0,})
        # formbase_0.update_default_vertex_attributes({'z': 0.0,})
        # gkey_key_base = formbase_0.gkey_key()

        # points_base = []
        dists = []
        i_s.append(i)

        # for key in formbase.vertices():
        #     points_base.append(formbase.vertex_coordinates(key))

        # for key in form.vertices_where({'is_fixed': False}):
        #     point = form.vertex_coordinates(key)
        #     try:
        #         gkey = geometric_key(point[:2] + [0])
        #         key_base = gkey_key_base[gkey]
        #         projected = formbase.vertex_coordinates(key_base)
        #         dist = distance_point_point(point, closest)
        #         dists.append(dist)
        #         form.set_vertex_attribute(key, name = 'dist', value = dist)
        #     except:
                
        #         dist, closest, _ = closest_point_in_cloud(point,points_base)
        #         form.set_vertex_attribute(key, name = 'dist', value = dist)
        #         dists.append(dist)

        for key in form.vertices():
            if form.get_vertex_attribute(key, 'is_fixed') == False:
                dists.append(form.get_vertex_attribute(key,'dist'))

        print('Form Optimised: {0}'.format(i))
        dist_max = max(dists)
        average = sum(dists)/len(dists)
        print('Maximum Distance is: {0}'.format(dist_max))
        print('Average Distance is: {0}'.format(average))
        averages.append(average)
        maxs.append(dist_max)

        plotter = MeshPlotter(form, figsize=(12, 8), tight=True)

        plotter.draw_vertices(
            keys=list(form.vertices_where({'is_external': False})),
            facecolor={key: i_to_white((attr['dist'] - 0) / (dist_max - 0)) for key, attr in form.vertices_where({'is_external': False}, True)},
            radius={key: 2 * attr['dist'] for key, attr in form.vertices_where({'is_external': False}, True)}
        )

        plotter.draw_edges(
            keys=list(form.edges_where({'is_edge': True})),
            color={key: '#00ff00' for key in form.edges_where({'is_external': True})},
            width={key: 2.0 for key in form.edges_where({'is_external': True})}
        )

        plotter.draw_faces(keys=list(form.faces_where({'is_loaded': True})))
        plotter.show()
        # plotter.save('/Users/mricardo/compas_dev/me/discretize/dists/img/0'+str(j)+'_0'+str(i)+'_dist.jpg')

    import matplotlib
    import matplotlib.pyplot as plt

    # Data for plotting
    print(i_s)
    print(averages)
    print(maxs)
    fig, ax = plt.subplots()
    ax.plot(i_s, averages)

    ax.set(xlabel='Discretization Level', ylabel='Average Distance (m)',
        title='Average Distance for Pattern 02')
    ax.grid()
    ax.set_xlim([0,10])
    ax.set_ylim([0,0.5])
    # fig.savefig("/Users/mricardo/compas_dev/me/discretize/dists/img/02_Average.png")
    plt.show()