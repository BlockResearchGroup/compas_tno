"""
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
    find_sym_axis_in_rect_patterns
    build_symmetry_matrix
    build_symmetry_transformation
    build_vertex_symmetry_transformation
    build_symmetry_matrix_supports

Constraints
====================

.. autosummary::
    :toctree: generated/

    check_envelope_constraints
    distance_target
    rectangular_smoothing_constraints
    assign_cracks
    rollers_on_openings
    set_b_constraint
    set_rollers_constraint

Loads
====================

.. autosummary::
    :toctree: generated/

    apply_selfweight_from_shape
    apply_selfweight_from_pattern
    apply_horizontal_multiplier
    apply_fill_load
    apply_selfweight_from_shape_proxy
    apply_selfweight_from_thrust


Envelopes
=========

.. autosummary::
    :toctree: generated/

    apply_envelope_from_shape
    apply_envelope_on_xy
    apply_envelope_on_xy_from_base
    apply_bounds_on_q
    project_mesh_to_middle
    modify_shapedata_with_spr_angle
    apply_envelope_from_shape_proxy
    apply_bounds_tub_tlb
    apply_bounds_reactions

Interpolation
====================

.. autosummary::
    :toctree: generated/

    interpolate_from_pointcloud
    get_shape_ub
    get_shape_ub_pattern
    get_shape_ub_fill
    get_shape_lb
    get_shape_lb_pattern
    get_shape_middle
    delaunay_mesh_from_points
    mesh_from_pointcloud
    create_mesh_from_topology_and_pointcloud
    create_mesh_from_topology_and_basemesh

Stiffness
====================

.. autosummary::
    :toctree: generated/

    compute_form_initial_lengths
    compute_edge_stiffness
    compute_average_edge_stiffness

Blocks
======

.. autosummary::
    :toctree: generated/

    blocks_from_dual
    extended_dual

Form Modifications
====================

.. autosummary::
    :toctree: generated/

    split_intersection_lines
    form_add_lines_support
    form_parabolic_slide
    move_pattern_to_origin
    fix_mesh_corners
    fix_mesh_boundary
    slide_diagram
    mesh_remove_two_valent_nodes
    slide_pattern_inwards
    displacement_map_paraboloid
    displacement_map_4parabolas
    displacement_map_parabola
    move_pattern_inwards
    move_pattern_outwards
    shuffle_diagram

Data Analysis
====================

.. autosummary::
    :toctree: generated/

    diagram_of_thrust
    diagram_of_multiple_thrust
    diagram_of_thrust_load_mult
    surface_GSF_load_mult
    save_csv_row
    open_csv_row
    interpolate_min_thk
    filter_min_thk
    lookup_folder
    save_pointcloud


"""
from __future__ import absolute_import

from .symmetry import (
    apply_radial_symmetry,
    apply_symmetry_from_axis,
    find_sym_axis_in_rect_patterns,
    build_symmetry_matrix,
    build_symmetry_transformation,
    build_vertex_symmetry_transformation,
    build_symmetry_matrix_supports,
)

from .constraints import (
    check_envelope_constraints,
    distance_target,
    rectangular_smoothing_constraints,
    assign_cracks,
    rollers_on_openings,
    set_b_constraint,
    set_rollers_constraint
)

from .functions import (
    paraboloid,
    dome
)

from .loads import (
    apply_selfweight_from_shape,
    apply_selfweight_from_pattern,
    apply_selfweight_from_thrust,
    apply_horizontal_multiplier,
    apply_fill_load,
    apply_selfweight_from_shape_proxy

)

from .envelopes import (
    apply_envelope_from_shape,
    apply_envelope_on_xy,
    apply_envelope_on_xy_from_base,
    apply_bounds_on_q,
    project_mesh_to_middle,
    modify_shapedata_with_spr_angle,
    apply_envelope_from_shape_proxy,
    apply_bounds_tub_tlb,
    apply_bounds_reactions
)

from .interpolation import (
    interpolate_from_pointcloud,
    get_shape_ub,
    get_shape_ub_pattern,
    get_shape_ub_fill,
    get_shape_lb,
    get_shape_lb_pattern,
    get_shape_middle,
    delaunay_mesh_from_points,
    mesh_from_pointcloud,
    create_mesh_from_topology_and_pointcloud,
    create_mesh_from_topology_and_basemesh
)

from .stiffness import (
    compute_form_initial_lengths,
    compute_edge_stiffness,
    compute_average_edge_stiffness
)

from .json import (
    update_json
)

from .blocks import (
    blocks_from_dual,
    extended_dual
)

from .form import (
    split_intersection_lines,
    form_add_lines_support,
    form_parabolic_slide,
    move_pattern_to_origin,
    fix_mesh_corners,
    fix_mesh_boundary,
    slide_diagram,
    mesh_remove_two_valent_nodes,
    store_inds,
    slide_pattern_inwards,
    displacement_map_paraboloid,
    displacement_map_4parabolas,
    displacement_map_parabola,
    move_pattern_inwards,
    move_pattern_outwards,
    shuffle_diagram
)

from .data_analysis import (
    diagram_of_thrust,
    diagram_of_multiple_thrust,
    diagram_of_thrust_load_mult,
    surface_GSF_load_mult,
    save_csv_row,
    open_csv_row,
    interpolate_min_thk,
    filter_min_thk,
    lookup_folder,
    save_pointcloud
)


__all__ = [

    'apply_radial_symmetry',
    'apply_symmetry_from_axis',
    'find_sym_axis_in_rect_patterns',
    'build_symmetry_matrix',
    'build_symmetry_transformation',
    'build_vertex_symmetry_transformation',
    'build_symmetry_matrix_supports',

    'check_envelope_constraints',
    'distance_target',
    'rectangular_smoothing_constraints',
    'assign_cracks',
    'rollers_on_openings',
    'set_b_constraint',
    'set_rollers_constraint',

    'paraboloid',
    'dome',

    'apply_selfweight_from_shape',
    'apply_selfweight_from_pattern',
    'apply_selfweight_from_thrust',
    'apply_horizontal_multiplier',
    'apply_fill_load',

    'apply_selfweight_from_shape_proxy',

    'apply_envelope_from_shape',
    'apply_envelope_on_xy',
    'apply_envelope_on_xy_from_base',
    'apply_bounds_on_q',
    'project_mesh_to_middle',
    'modify_shapedata_with_spr_angle',
    'apply_envelope_from_shape_proxy',
    'apply_bounds_tub_tlb',
    'apply_bounds_reactions',

    'interpolate_from_pointcloud',
    'get_shape_ub',
    'get_shape_ub_pattern',
    'get_shape_ub_fill',
    'get_shape_lb',
    'get_shape_lb_pattern',
    'get_shape_middle',
    'delaunay_mesh_from_points',
    'mesh_from_pointcloud',
    'create_mesh_from_topology_and_pointcloud',
    'create_mesh_from_topology_and_basemesh',

    'compute_form_initial_lengths',
    'compute_edge_stiffness',
    'compute_average_edge_stiffness',

    'update_json',

    'blocks_from_dual',
    'extended_dual',

    'split_intersection_lines',
    'form_add_lines_support',
    'form_parabolic_slide',
    'move_pattern_to_origin',
    'fix_mesh_corners',
    'fix_mesh_boundary',
    'slide_diagram',
    'mesh_remove_two_valent_nodes',
    'store_inds',
    'slide_pattern_inwards',
    'displacement_map_paraboloid',
    'displacement_map_4parabolas',
    'displacement_map_parabola',
    'move_pattern_inwards',
    'move_pattern_outwards',
    'shuffle_diagram',

    'diagram_of_thrust',
    'diagram_of_multiple_thrust',
    'diagram_of_thrust_load_mult',
    'surface_GSF_load_mult',
    'save_csv_row',
    'open_csv_row',
    'interpolate_min_thk',
    'filter_min_thk',
    'lookup_folder',
    'save_pointcloud',
]
