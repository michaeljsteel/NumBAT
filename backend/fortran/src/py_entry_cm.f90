! Let python know the compiled version of the code
! To check git changes that have been pulled but not built
subroutine get_fortran_compiled_nb_version(ver_str)
   use numbatmod

   character(len=EMSG_LENGTH), intent(out) :: ver_str

   write(ver_str, '(A)') NUMBAT_VERSION_STR_MMM

   end subroutine

subroutine calc_em_modes(n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    !  inputs
   E_H_field, bdy_cdn, itermax, arp_tol, debug, &
   mesh_file, n_msh_pts, n_msh_elts, n_elt_mats, v_refindex_n, shortrun, & !  inputs
   v_evals_beta, femsol_em, poln_fracs, &
   m_elnd_to_mshpt, v_elt_material, v_mshpt_physindex, v_mshpt_xy, &
   ls_material, errco, emsg)

   use numbatmod
   use calc_em_impl


   integer(8), intent(in) :: n_modes
   double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
   complex(8), intent(in) :: shift_ksqr

   integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
   double precision, intent(in) :: arp_tol
   character(len=*), intent(in) :: mesh_file
   integer(8), intent(in) :: n_msh_pts,  n_msh_elts, n_elt_mats

   complex(8), intent(in) ::  v_refindex_n(n_elt_mats)
   integer(8) :: shortrun

   complex(8), intent(out) :: v_evals_beta(n_modes)
   complex(8), intent(out) :: femsol_em(3,P2_NODES_PER_EL+7,n_modes,n_msh_elts)

   complex(8), intent(out) :: poln_fracs(4,n_modes)
   integer(8), intent(out) :: m_elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   integer(8), intent(out) :: v_elt_material(n_msh_elts), v_mshpt_physindex(n_msh_pts)
   double precision, intent(out) :: v_mshpt_xy(2,n_msh_pts)
   complex(8), intent(out) :: ls_material(1,P2_NODES_PER_EL+7,n_msh_elts)

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg
   type(NBError) nberr


     call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

     !call nberr%reset()

     ! arp_tol = 1.0d-12 ! TODO: ARPACK_ stopping precision,  connect  to user switch

     call calc_em_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &
     E_H_field, bdy_cdn, itermax, arp_tol, debug, mesh_file,&
     n_msh_pts, n_msh_elts, n_elt_mats, v_refindex_n, shortrun, &
     v_evals_beta, femsol_em, poln_fracs, m_elnd_to_mshpt, &
     v_elt_material, v_mshpt_physindex, v_mshpt_xy, ls_material, nberr)

    ! call nberr%to_py(errco, emsg)

     call finalise_pyerror(nberr, errco, emsg)

   end subroutine

subroutine calc_ac_modes(n_modes, q_ac, dimscale_in_m, shift_nu, &
   bdy_cdn, itermax, arp_tol,  &
   symmetry_flag, n_elt_mats, c_tensor, rho, build_acmesh_from_emmesh, &
   mesh_file, n_msh_pts, n_msh_elts, &
   v_mshpt_physindex, &
   m_elnd_to_mshpt, v_elt_material, v_mshpt_xy, &
   v_evals_nu, femsol_ac, poln_fracs, errco, emsg)

   use numbatmod
   use calc_ac_impl

   integer(8), intent(in) :: n_modes

   complex(8), intent(in) :: q_ac
   double precision, intent(in) :: dimscale_in_m
   complex(8), intent(in) :: shift_nu
   integer(8), intent(in) :: bdy_cdn, itermax
   double precision, intent(in) :: arp_tol
   integer(8), intent(in) :: symmetry_flag, build_acmesh_from_emmesh
   integer(8), intent(in) :: n_elt_mats

   complex(8), intent(in) :: c_tensor(6,6,n_elt_mats)
   complex(8), intent(in) :: rho(n_elt_mats)

   character(len=FNAME_LENGTH), intent(in)  :: mesh_file
   integer(8), intent(in) :: n_msh_pts, n_msh_elts

   integer(8), intent(in) :: v_mshpt_physindex(n_msh_pts)

   integer(8) :: v_elt_material(n_msh_elts)
   integer(8) :: m_elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)

   double precision ::  v_mshpt_xy(2,n_msh_pts)

   complex(8), intent(out) :: v_evals_nu(n_modes)
   complex(8), intent(out) :: femsol_ac(3,P2_NODES_PER_EL,n_modes,n_msh_elts)

   complex(8), intent(out) :: poln_fracs(4,n_modes)

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   type(NBError) nberr

   !call nberr%reset()
     call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls


   !f2py intent(in) m_elnd_to_mshpt, v_elt_material, v_mshpt_xy
   !f2py intent(out) m_elnd_to_mshpt, v_elt_material, v_mshpt_xy

   call calc_ac_modes_impl(n_modes, q_ac, dimscale_in_m, shift_nu, &
      bdy_cdn, itermax, arp_tol, &
      symmetry_flag, c_tensor, rho, build_acmesh_from_emmesh, &
      mesh_file, n_msh_pts, n_msh_elts, n_elt_mats,  &
      m_elnd_to_mshpt, v_elt_material, v_mshpt_physindex, v_mshpt_xy, &
      v_evals_nu, femsol_ac, poln_fracs, nberr)

   !call nberr%to_py(errco, emsg)
   call finalise_pyerror(nberr, errco, emsg)

end subroutine

subroutine conv_gmsh(geo_name, assertions_on, errco, emsg)

   use numbatmod

   !TODO: f2py doesn't seem to be able to handle this length being set by the #include and preprocessor
   integer(8) :: assertions_on
   character(len=*), intent(in) :: geo_name

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   type(NBError) nberr

   !call nberr%reset()

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call conv_gmsh_impl(geo_name, assertions_on, nberr)

   !call nberr%to_py(errco, emsg)

   call finalise_pyerror(nberr, errco, emsg)

end
