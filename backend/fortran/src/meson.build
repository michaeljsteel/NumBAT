fortmodc_h = files('myhead.h')

nb_module_fsrcs = files(
    'alloc_module.f90',
    'numbatmod_module.f90',
    'stopwatch_module.f90',
     'conv_gmsh_interfaces.f90',

    #'nbinterfaces_module.f90',
)

nb_pyinterface_fsrcs = files(
    'f2py_entry_functions.f90',

    'array_material_em.f90',
    'array_material_ac.f90',


    'py_entry_cm.f90',
)


nb_impl_fsrcs = files(
    'array_size.f90',

    'ac_alpha_analytic.f90',
    'ac_alpha_quadrature.f90',

    'ac_mode_energy_analytic.f90',
    'ac_mode_energy_quadrature.f90',

    'ac_mode_power_analytic.f90',
    'ac_mode_power_quadrature.f90',

    'em_mode_energy_int.f90',
    'em_mode_act_energy_quadrature.f90',

    'em_mode_power_sz_quadrature.f90',
    'em_mode_power_sz_analytic.f90',

    'h_mode_field_ez.f90',

    'moving_boundary.f90',
    'photoelastic_int_common.f90',


  #'quicksort.f90',
  'construct_solution_fields_ac.f90',
  'construct_solution_fields_em.f90',
  #'array_sol.f90',
    'ac_mode_power_analytic_matel_powerflow.f90',

    'build_fem_ops_ac.f90',
    'build_fem_ops_em.f90',
    'basis_functions.f90',
    'bound_cond_ac.f90',
    'construct_fem_tables.f90',
    'conv_gmsh_subroutines_mod.f90',
    'csr_length_AC.f90',
    'conv_gmsh_impl.f90',

    #'csr_max_length_ac.f90',
    'get_H_field_p3.f90',
    'get_coord_p3.f90',
    'gmsh_post_process.f90',
    'gmsh_post_process_ac.f90',
    'h_mode_field_ez.f90',
    'interp_nod_2d.f90',
    'is_curved_elem_tri.f90',
    'jacobian_p1_2d.f90',
    'jacobian_p2_2d.f90',
    'mat_el.f90',
    'mat_el_powerflow_old.f90',
    'mat_el_energy_rho.f90',
    'mat_el_v2.f90',
    'make_elt_femops_ac_no_sym.f90',
    'mat_p2_p2.f90',
    'mat_p2_p2_p2.f90',
    'mat_p2_p2_p2_1d.f90',
    'mat_p2_p2_p2x.f90',
    'mat_p2_p2_p2y.f90',
    'mat_p2_p2x.f90',
    'mat_p2_p2y.f90',
    'mat_p2x_p2x.f90',
    'mat_p2x_p2y.f90',
    'mat_p2y_p2y.f90',
    'meshprops.f90',
    'mode_energy.f90',
    'normalisation.f90',
    'orthogonal.f90',
    'periodic_N_E_F.f90',
    'periodic_cond.f90',
    'periodic_lattice_vec.f90',
    'periodic_node.f90',
    'phi1_2d_mat.f90',
    'phi2_2d_mat.f90',
    'phi2_2d_mat_J.f90',
    'phi2_grad.f90',
    'phi3_2d_mat.f90',
    'phi3_grad_p2.f90',
    #'phi_vector_map.f90',
   # 'phi_vector_basis.f90',
#'photoelastic_int_curvilinear_impl.f90',
#    'photoelastic_int_linear_impl.f90',
    'photoelastic_int_common_impl.f90',

    'calc_modes_em_impl.f90',
    'calc_modes_ac_impl.f90',
    'periodicbcs_class.f90',
    'quad_triangle.f90',
    'quad_integrator.f90',
    'sparse_csc_common.f90',
    'sparse_csc_em.f90',
    'sparse_csc_ac.f90',
    'sort_csc.f90',
    'sort_int.f90',
    #'sort_n.f90',
    'utils.f90',
    'solve_arpack_problem.f90',
    #'valpr_umf_64_ac.f90',
    'find_eigvals_order.f90',
    #'find_eigvals_order_AC.f90',
    'z_mxv_csc_64.f90',

    #'array_size.f90',
    #'bound_cond_em.f90',
    #'build_mesh_tables.f90',
    #'conv_gmsh.f90',
    #'conv_gmsh_symmetry.f90',
    #'csr_make_arrays.f90',
    #'csr_loose.f90',
    #'label_edges.f90',
    #'label_faces.f90',
    #'label_nodes_P3.f90',
    #'py_calc_modes_impl.f90',
    #'type_node_edge_face.f90',
)


nb_pyf = custom_target(
    'f2py_make_pyf',
    input: nb_pyinterface_fsrcs,
    output: ['nb_fortran.pyf'],
    command: [
        comp_f2py,
        '-m', 'nb_fortran',
        '-h', 'nb_fortran.pyf',
        '--overwrite-signature',
        '--lower', '@INPUT@',
    ],
)


nb_fortmodc = custom_target(
    'f2py_make_cmod',
    input: [nb_pyinterface_fsrcs],
    output: ['nb_fortranmodule.c',
    'nb_fortran-f2pywrappers.f'],
    command: [
        comp_f2py,
        '-m', 'nb_fortran',
        '--build-dir', 'src',
        '--lower',
        '-include"@0@"'.format(
         meson.current_source_dir() / 'myhead.h'),
        '@INPUT@',
    ],
)

