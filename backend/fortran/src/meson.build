nb_mod_fsrcs = files([ 'numbatmod_module.f90', 'nbinterfaces_module.f90', 'alloc_module.f90', 'stopwatch_module.f90', ])

nb_pyint_fsrcs = files(['pyinterfaces.f90', 'conv_gmsh.f90', 'em_mode_energy_int_v2_ez.f',
                  'em_mode_energy_int_ez.f', 'em_mode_e_energy_int.f',
                  'h_mode_field_ez.f90', 'ac_mode_power_int.f', 'ac_mode_power_int_v2.f',
                  'ac_mode_power_int_v4.f', 'ac_mode_elastic_energy_int.f',
                  'ac_mode_elastic_energy_int_v4.f', 'ac_alpha_int.f',
                  'ac_alpha_int_v2.f', 'array_size.f90', 'array_material_ac.f',
                  'photoelastic_int.f90', 'photoelastic_int_v2.f90',
                  'moving_boundary.f90', 'py_entry_cm.f90'])

#nb_fsrcs_tt = ['src/array_size.f90']
nb_fsrcs = files( ['array_size.f90', 'array_sol_ac_v2.f90', 'array_sol.f90',
            'asmbly_ac.f90', 'asmbly.f90', 'basis_ls.f90', 'basis_vec.f90',
            'bound_cond_ac.f90', 'bound_cond.f90', 'build_mesh_tables.f90',
            'construct_fem_tables.f90', 'conv_gmsh.f90', 'conv_gmsh_subroutines_mod.f90',
            'conv_gmsh_symmetry.f90', 'csr_make_arrays.f90', 'csr_max_length_ac.f90',
            'csr_max_length.f90', 'get_coord_p3.f90', 'get_H_field_p3.f90',
            'gmsh_post_process_ac.f90', 'gmsh_post_process.f90', 'h_mode_field_ez.f90',
            'interp_nod_2d.f90', 'is_curved_elem_tri.f90', 'jacobian_p1_2d.f90',
            'jacobian_p2_2d.f90', 'label_edges.f90', 'label_faces.f90',
            'label_nodes_P3.f90', 'lattice_vec.f90', 'mat_el_energy.f90',
            'mat_el_energy_rho.f90', 'mat_el.f90', 'mat_el_v2.f90', 'mat_el_v3.f90',
            'mat_p2_p2.f90', 'mat_p2_p2_p2_1d.f90', 'mat_p2_p2_p2.f90', 'mat_p2_p2_p2x.f90',
            'mat_p2_p2_p2y.f90', 'mat_p2_p2x.f90', 'mat_p2_p2y.f90', 'mat_p2x_p2x.f90',
            'mat_p2x_p2y.f90', 'mat_p2y_p2y.f90', 'meshprops.f90', 'mode_energy.f90',
            'moving_boundary.f90', 'normalisation.f90', 'orthogonal.f90', 'periodic_cond.f90',
            'periodic_N_E_F.f90', 'periodic_node.f90', 'phi1_2d_mat.f90', 'phi2_2d_mat.f90',
            'phi2_2d_mat_J.f90', 'phi2_grad.f90', 'phi3_2d_mat.f90', 'phi3_grad_p2.f90',
            'photoelastic_int.f90', 'photoelastic_int_v2.f90', 'py_calc_modes_ac.f90',
            'py_calc_modes.f90', 'py_calc_modes_impl.f90', 'py_entry_cm.f90',
            'pyinterfaces.f90', 'quad_triangle.f90', 'sort_csr.f90', 'sort_int.f90',
            'sort_n.f90', 'type_node_edge_face.f90', 'utils.f90', 'valpr_umf_64_ac.f90',
            'valpr_umf_64.f90', 'z_indexx_AC.f90', 'z_indexx.f90', 'z_mxv_csc_64.f90'
            ])


comp_f2py = find_program('f2py')

nb_pyf1 = custom_target('f2py_out_stage1',
                        input : nb_pyint_fsrcs,
                        output : ['nb_fortran.pyf'],
                        command : [comp_f2py,
                                  '-m' ,'nb_fortran','-h', 'nb_fortran.pyf',
                                  '--overwrite-signature', '--lower',
                                  '@INPUT@'
                                  ])
nb_pyf2 = custom_target('f2py_out_stage2',
                        input: nb_pyint_fsrcs,
                        output: [
                        'nb_fortranmodule.c'],
                        command: [comp_f2py,
                                  '-m', 'nb_fortran', '--build-dir', 'f2pybuild',
                                  '--lower', '-include \"../src/myhead.h\" ',
                                  '@INPUT@'
                                  ])