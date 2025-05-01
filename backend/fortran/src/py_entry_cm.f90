
subroutine calc_em_modes(n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    !  inputs
   E_H_field, bdy_cdn, itermax, debug, &
   mesh_file, n_msh_pts, n_msh_el, n_elt_mats, v_refindex_n, shortrun, & !  inputs
   v_evals_beta, femsol_em, poln_fracs, &
   elnd_to_mshpt, v_el_material, v_nd_physindex, v_nd_xy, &
   ls_material, errco, emsg)

   use numbatmod
   use calc_em_impl


   integer(8), intent(in) :: n_modes
   double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
   complex(8), intent(in) :: shift_ksqr

   integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
   character(len=*), intent(in) :: mesh_file
   integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_elt_mats

   complex(8), intent(in) ::  v_refindex_n(n_elt_mats)
   integer(8) :: shortrun

   complex(8), intent(out) :: v_evals_beta(n_modes)
   complex(8), intent(out) :: femsol_em(3,P2_NODES_PER_EL+7,n_modes,n_msh_el)

   complex(8), intent(out) :: poln_fracs(4,n_modes)
   integer(8), intent(out) :: elnd_to_mshpt(P2_NODES_PER_EL, n_msh_el)
   integer(8), intent(out) :: v_el_material(n_msh_el), v_nd_physindex(n_msh_pts)
   double precision, intent(out) :: v_nd_xy(2,n_msh_pts)
   complex(8), intent(out) :: ls_material(1,P2_NODES_PER_EL+7,n_msh_el)
   double precision arp_tol

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg
   type(NBError) nberr

   call nberr%reset()

 arp_tol = 1.0d-12 ! TODO: ARPACK_ stopping precision,  connect  to user switch

   call calc_em_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &
      E_H_field, bdy_cdn, itermax, arp_tol, debug, mesh_file,&
      n_msh_pts, n_msh_el, n_elt_mats, v_refindex_n, shortrun, &
      v_evals_beta, femsol_em, poln_fracs, elnd_to_mshpt, &
      v_el_material, v_nd_physindex, v_nd_xy, ls_material, nberr)

   call nberr%to_py(errco, emsg)

end subroutine

subroutine calc_ac_modes(n_modes, q_ac, dimscale_in_m, shift_nu, &
   bdy_cdn, itermax, arp_tol, debug, &
   symmetry_flag, n_elt_mats, c_tensor, rho, build_acmesh_from_emmesh, &
   mesh_file, n_msh_pts, n_msh_el, &
   v_nd_physindex, &
   elnd_to_mshpt, v_el_material, v_nd_xy, &
   v_eigs_nu, femsol_ac, poln_fracs, errco, emsg)

   use numbatmod
   use calc_ac_impl

   integer(8), intent(in) :: n_modes

   complex(8), intent(in) :: q_ac
   double precision, intent(in) :: dimscale_in_m
   complex(8), intent(in) :: shift_nu
   integer(8), intent(in) :: bdy_cdn, itermax, debug
   double precision, intent(in) :: arp_tol
   integer(8), intent(in) :: symmetry_flag, build_acmesh_from_emmesh
   integer(8), intent(in) :: n_elt_mats

   complex(8), intent(in) :: c_tensor(6,6,n_elt_mats)
   complex(8), intent(in) :: rho(n_elt_mats)

   character(len=FNAME_LENGTH), intent(in)  :: mesh_file
   integer(8), intent(in) :: n_msh_pts, n_msh_el

   integer(8), intent(in) :: v_nd_physindex(n_msh_pts)

   integer(8) :: v_el_material(n_msh_el)
   integer(8) :: elnd_to_mshpt(P2_NODES_PER_EL, n_msh_el)

   double precision ::  v_nd_xy(2,n_msh_pts)

   complex(8), intent(out) :: v_eigs_nu(n_modes)
   complex(8), intent(out) :: femsol_ac(3,P2_NODES_PER_EL,n_modes,n_msh_el)

   complex(8), intent(out) :: poln_fracs(4,n_modes)

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   type(NBError) nberr

   call nberr%reset()

   !f2py intent(in) elnd_to_mshpt, v_el_material, v_nd_xy
   !f2py intent(out) elnd_to_mshpt, v_el_material, v_nd_xy

   call calc_ac_modes_impl(n_modes, q_ac, dimscale_in_m, shift_nu, &
      bdy_cdn, itermax, arp_tol, debug,  &
      symmetry_flag, c_tensor, rho, build_acmesh_from_emmesh, &
      mesh_file, n_msh_pts, n_msh_el, n_elt_mats,  &
      elnd_to_mshpt, v_el_material, v_nd_physindex, v_nd_xy, &
      v_eigs_nu, femsol_ac, poln_fracs, nberr)

   call nberr%to_py(errco, emsg)

end subroutine
