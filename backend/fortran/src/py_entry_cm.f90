
subroutine calc_em_modes( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    !  inputs
        E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_elt_mats, &
        v_refindex_n, & !  inputs
        v_eigs_beta, sol1, mode_pol, elnd_to_mesh, type_el, type_nod, &
        v_nd_xy, ls_material, errco, emsg)

    use numbatmod
    use calc_em_impl


    integer(8), parameter  :: d_nodes_per_el = 6
    integer(8), intent(in) :: n_modes
    double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
    complex(8), intent(in) :: shift_ksqr

    integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
    character(len=*), intent(in) :: mesh_file
    integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_elt_mats

    complex(8), intent(in) ::  v_refindex_n(n_elt_mats)

    complex(8), intent(out) :: v_eigs_beta(n_modes)
    complex(8), intent(out) :: sol1(3,d_nodes_per_el+7,n_modes,n_msh_el)

    complex(8), intent(out) :: mode_pol(4,n_modes)
    integer(8), intent(out) :: elnd_to_mesh(d_nodes_per_el, n_msh_el)
    integer(8), intent(out) :: type_el(n_msh_el), type_nod(n_msh_pts)
    double precision, intent(out) :: v_nd_xy(2,n_msh_pts)
    complex(8), intent(out) :: ls_material(1,d_nodes_per_el+7,n_msh_el)

    integer(8),  intent(out) :: errco
    character(len=EMSG_LENGTH), intent(out) :: emsg

    call calc_em_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &
        E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_elt_mats, v_refindex_n, &
        v_eigs_beta, sol1, mode_pol, elnd_to_mesh, type_el, type_nod, v_nd_xy, ls_material, errco, emsg)

end subroutine

subroutine calc_ac_modes(n_modes, q_ac, dimscale_in_m, shift_nu, &
        i_bnd_cdns, itermax, tol, debug, show_mem_est, &
        symmetry_flag, n_elt_mats, c_tensor, rho, supplied_geo_flag, &
        mesh_file, n_msh_pts, n_msh_el, &
        type_nod, &
        elnd_to_mesh, type_el, v_nd_xy, &
        v_eigs_nu, sol1, mode_pol, errco, emsg)

    use numbatmod
    use calc_ac_impl

    integer(8),  parameter :: d_nodes_per_el = 6

    integer(8), intent(in) :: n_modes

    complex(8), intent(in) :: q_ac
    double precision, intent(in) :: dimscale_in_m
    complex(8), intent(in) :: shift_nu
    integer(8), intent(in) :: i_bnd_cdns, itermax, debug, show_mem_est
    double precision, intent(in) :: tol
    integer(8), intent(in) :: symmetry_flag, supplied_geo_flag
    integer(8), intent(in) :: n_elt_mats

    complex(8), intent(in) :: c_tensor(6,6,n_elt_mats)
    complex(8), intent(in) :: rho(n_elt_mats)

    character(len=FNAME_LENGTH), intent(in)  :: mesh_file
    integer(8), intent(in) :: n_msh_pts, n_msh_el

    integer(8), intent(in) :: type_nod(n_msh_pts)

    integer(8) :: type_el(n_msh_el)
    integer(8) :: elnd_to_mesh(d_nodes_per_el, n_msh_el)

    double precision ::  v_nd_xy(2,n_msh_pts)

    complex(8), intent(out) :: v_eigs_nu(n_modes)
    complex(8), intent(out) :: sol1(3,d_nodes_per_el,n_modes,n_msh_el)

    complex(8), intent(out) :: mode_pol(4,n_modes)

    integer(8),  intent(out) :: errco
    character(len=EMSG_LENGTH), intent(out) :: emsg

    !f2py intent(in) elnd_to_mesh, type_el, v_nd_xy
    !f2py intent(out) elnd_to_mesh, type_el, v_nd_xy

    call calc_ac_modes_impl(n_modes, q_ac, dimscale_in_m, shift_nu, &
        i_bnd_cdns, itermax, tol, debug, show_mem_est, &
        symmetry_flag, n_elt_mats, c_tensor, rho, supplied_geo_flag, &
        mesh_file, n_msh_pts, n_msh_el, &
        type_nod, &
        elnd_to_mesh, type_el, v_nd_xy, &
        v_eigs_nu, sol1, mode_pol, errco, emsg)

end subroutine
