
subroutine calc_EM_modes( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    ! inputs
        E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_typ_el, &
        v_refindex_n, & ! inputs
        v_eigs_beta, sol1, mode_pol, table_nod, type_el, type_nod, &
        mesh_xy, ls_material, errco, emsg)  

    use numbatmod
    use calc_em_impl


    integer(8), parameter  :: d_nodes_per_el = 6
    integer(8), intent(in) :: n_modes
    double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
    complex(8), intent(in) :: shift_ksqr

    integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
    character(len=*), intent(in) :: mesh_file
    integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_typ_el

    complex(8), intent(in) ::  v_refindex_n(n_typ_el)
    !complex(8), intent(in) ::  v_refindex_n(*)

    complex(8), intent(out) :: v_eigs_beta(n_modes)
    complex(8), intent(out) :: sol1(3,d_nodes_per_el+7,n_modes,n_msh_el)

    complex(8), intent(out) :: mode_pol(4,n_modes)
    integer(8), intent(out) :: table_nod(d_nodes_per_el, n_msh_el)
    integer(8), intent(out) :: type_el(n_msh_el), type_nod(n_msh_pts)
    double precision, intent(out) :: mesh_xy(2,n_msh_pts)
    complex(8), intent(out) :: ls_material(1,d_nodes_per_el+7,n_msh_el)

    integer, intent(out) :: errco
    character(len=EMSG_LENGTH), intent(out) :: emsg

    write(*,*) 'calc_em_ 1'
    write(*,*) v_refindex_n
    write(*,*) 'calc_em_ 2'
    write(*,*) 'calc_em_ 3'
    call calc_EM_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    
        E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_typ_el, v_refindex_n, & 
        v_eigs_beta, sol1, mode_pol, table_nod, type_el, type_nod, mesh_xy, ls_material, errco, emsg) 

end subroutine

subroutine calc_ac_modes(n_modes, q_ac, dimscale_in_m, shift_nu, &
        i_bnd_cdns, itermax, tol, debug, show_mem_est, &
        symmetry_flag, n_typ_el, c_tensor, rho, supplied_geo_flag, &
        mesh_file, n_msh_pts, n_msh_el, &
        type_nod, &
        table_nod, type_el, mesh_xy, &
        v_eigs_nu, sol1, mode_pol, errco, emsg)

    use numbatmod
    use calc_ac_impl

    integer, parameter :: d_nodes_per_el = 6

    integer(8), intent(in) :: n_modes

    complex(8), intent(in) :: q_ac
    double precision, intent(in) :: dimscale_in_m
    complex(8), intent(in) :: shift_nu
    integer(8), intent(in) :: i_bnd_cdns, itermax, debug, show_mem_est
    double precision, intent(in) :: tol
    integer(8), intent(in) :: symmetry_flag, supplied_geo_flag
    integer(8), intent(in) :: n_typ_el

    complex(8), intent(in) :: c_tensor(6,6,n_typ_el)
    complex(8), intent(in) :: rho(n_typ_el)

    character(len=FNAME_LENGTH), intent(in)  :: mesh_file
    integer(8), intent(in) :: n_msh_pts, n_msh_el

    integer(8), intent(in) :: type_nod(n_msh_pts)

    integer(8) :: type_el(n_msh_el)
    integer(8) :: table_nod(d_nodes_per_el, n_msh_el)

    double precision ::  mesh_xy(2,n_msh_pts)

    complex(8), intent(out) :: v_eigs_nu(n_modes)
    complex(8), intent(out) :: sol1(3,d_nodes_per_el,n_modes,n_msh_el)

    complex(8), intent(out) :: mode_pol(4,n_modes)

    integer, intent(out) :: errco
    character(len=EMSG_LENGTH), intent(out) :: emsg

    !f2py intent(in) table_nod, type_el, mesh_xy
    !f2py intent(out) table_nod, type_el, mesh_xy

    call calc_ac_modes_impl(n_modes, q_ac, dimscale_in_m, shift_nu, &
        i_bnd_cdns, itermax, tol, debug, show_mem_est, &
        symmetry_flag, n_typ_el, c_tensor, rho, supplied_geo_flag, &
        mesh_file, n_msh_pts, n_msh_el, &
        type_nod, &
        table_nod, type_el, mesh_xy, &
        v_eigs_nu, sol1, mode_pol, errco, emsg)

end subroutine
