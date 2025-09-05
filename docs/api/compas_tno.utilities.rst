********************************************************************************
compas_tno.utilities
********************************************************************************

.. currentmodule:: compas_tno.utilities

Symmetry
====================

.. autosummary::
    :toctree: generated/

    apply_radial_symmetry
    apply_symmetry_from_axis
    build_symmetry_matrix
    build_symmetry_matrix_supports
    build_symmetry_transformation
    build_vertex_symmetry_transformation
    find_sym_axis_in_rect_patterns

Constraints
====================

.. autosummary::
    :toctree: generated/

    assign_cracks
    check_envelope_constraints
    distance_target
    rectangular_smoothing_constraints
    rollers_on_openings
    set_b_constraint
    set_rollers_constraint


Interpolation
====================

.. autosummary::
    :toctree: generated/

    create_mesh_from_topology_and_basemesh
    create_mesh_from_topology_and_pointcloud
    interpolate_from_pointcloud
    mesh_from_pointcloud

Stiffness
====================

.. autosummary::
    :toctree: generated/

    compute_average_edge_stiffness
    compute_edge_stiffness
    compute_form_initial_lengths


Form Modifications
====================

.. autosummary::
    :toctree: generated/

    displacement_map_4parabolas
    displacement_map_parabola
    displacement_map_paraboloid
    fix_mesh_boundary
    fix_mesh_corners
    form_add_lines_support
    form_parabolic_slide
    mesh_remove_two_valent_nodes
    move_pattern_inwards
    move_pattern_outwards
    move_pattern_to_origin
    shuffle_diagram
    slide_diagram
    slide_pattern_inwards
    split_intersection_lines
    store_inds

