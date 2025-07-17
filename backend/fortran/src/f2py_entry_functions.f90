#include "numbat_decl.h"
! Calculate the v_alpha integral of an AC mode with itself using
! Direct integration

! \alpha = \Omega^2/Energy_aC \int  eta_ijkl d_i u_j^* d_k u_l

subroutine ac_alpha_analytic (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material,  &
   eta_ijkl, q_AC, Omega_AC, soln_ac_u, &
   v_ac_mode_energy, v_alpha_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC, v_ac_mode_energy(n_modes)
   complex(8) eta_ijkl(3,3,3,3,n_elt_mats)
   double precision, dimension(n_modes), intent(out) :: v_alpha_r
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, eta_ijkl, q_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC, v_ac_mode_energy

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(eta_ijkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!f2py depend(v_ac_mode_energy) n_modes

    type(NBError) nberr

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call ac_alpha_analytic_impl (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material,  &
   eta_ijkl, q_AC, Omega_AC, soln_ac_u, v_ac_mode_energy, v_alpha_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)


end subroutine ac_alpha_analytic


subroutine ac_alpha_quadrature (n_modes, n_msh_elts, n_msh_pts, &
   v_mshpt_xy, m_elnd_to_mshpt, n_elt_mats, v_elt_material, eta_ijkl, &
   q_AC, Omega_AC, soln_ac_u, &
   v_ac_mode_energy, v_alpha_r, errco, emsg)


   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes
   integer(8) n_msh_elts, n_msh_pts,  n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC, v_ac_mode_energy(n_modes)
   complex(8) eta_ijkl(3,3,3,3,n_elt_mats)

   double precision, dimension(n_modes), intent(out) :: v_alpha_r
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, eta_ijkl, q_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC, v_ac_mode_energy

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(eta_ijkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!f2py depend(v_ac_mode_energy) n_modes

!f2py intent(out) v_alpha


    type(NBError) nberr

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call ac_alpha_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, &
   v_mshpt_xy, m_elnd_to_mshpt, n_elt_mats, v_elt_material, eta_ijkl, &
   q_AC, Omega_AC, soln_ac_u, v_ac_mode_energy, v_alpha_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)


end subroutine ac_alpha_quadrature

subroutine ac_mode_energy_analytic (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) rho(n_elt_mats)

   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)

   double precision, dimension(n_modes), intent(out) :: v_energy_r

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, rho
!f2py intent(in) soln_ac_u, Omega_AC

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(Omega_AC) n_modes
!f2py depend(rho) n_elt_mats

    type(NBError) nberr

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call ac_mode_energy_analytic_impl (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine ac_mode_energy_analytic


subroutine AC_mode_energy_quadrature (n_modes, n_msh_elts, n_msh_pts, v_mshpt_xy, &
   m_elnd_to_mshpt, n_elt_mats, v_elt_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) Omega_AC(n_modes)
   double precision, dimension(n_modes), intent(out) :: v_energy_r

   complex(8) rho(n_elt_mats)

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   type(NBError) nberr


!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, rho
!f2py intent(in) soln_ac_u, debug, Omega_AC
!
!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(Omega_AC) n_modes
!f2py depend(rho) n_elt_mats

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call AC_mode_energy_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, v_mshpt_xy, &
   m_elnd_to_mshpt, n_elt_mats, v_elt_material,  &
   rho, Omega_AC, soln_ac_u, v_energy_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine AC_mode_energy_quadrature





subroutine ac_mode_power_analytic (n_modes, n_msh_elts, n_msh_pts,  &
   v_mshpt_xy, m_elnd_to_mshpt, &
   n_elt_mats, v_elt_material, stiffC_IJ_el, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz_r, errco, emsg)

   use numbatmod
   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, n_msh_elts, n_msh_pts
   double precision v_mshpt_xy(2,n_msh_pts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)

   integer(8) n_elt_mats
   integer(8) v_elt_material(n_msh_elts)

   complex(8) stiffC_IJ_el(6,6,n_elt_mats)

   complex(8) q_AC
   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   double precision, dimension(n_modes), intent(out) :: v_power_Sz_r
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   type(NBError) nberr


!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, stiffC_zjkl, q_AC
!f2py intent(in) soln_ac_u, Omega_AC
!
!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes

!f2py intent(out) v_power_Sz

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call ac_mode_power_analytic_impl (n_modes, n_msh_elts, n_msh_pts,  &
   v_mshpt_xy, m_elnd_to_mshpt, &
   n_elt_mats, v_elt_material, stiffC_IJ_el, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine ac_mode_power_analytic



subroutine ac_mode_power_quadrature (n_modes, n_msh_elts, n_msh_pts, &
   v_mshpt_xy, m_elnd_to_mshpt,  &
   n_elt_mats, v_elt_material, stiffC_zjkl, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz_r, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, n_msh_elts, n_msh_pts
   double precision v_mshpt_xy(2,n_msh_pts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   integer(8) n_elt_mats, v_elt_material(n_msh_elts)

   complex(8) stiffC_zjkl(3,3,3,n_elt_mats)
   complex(8) q_AC
   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   double precision, dimension(n_modes) :: v_power_Sz_r
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg




!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, stiffC_zjkl, q_AC
!f2py intent(in) soln_ac_u, Omega_AC

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes

!f2py intent(out) v_power_Sz_r

   type(NBError) nberr

 call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call ac_mode_power_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, &
   v_mshpt_xy, m_elnd_to_mshpt,  &
   n_elt_mats, v_elt_material, stiffC_zjkl, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz_r, nberr)

   call finalise_pyerror(nberr, errco, emsg)


end subroutine ac_mode_power_quadrature



subroutine EM_mode_energy_int (k_0, n_modes, n_msh_elts, n_msh_pts,&
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_energy, errco, emsg)
   !
   !     k_0 = 2 pi / lambda, where lambda in meters.
   !
   use numbatmod
   use class_TriangleIntegrators

   integer(8) n_modes, n_msh_elts, n_msh_pts
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_em_e(3,N_DOF_PER_EL,n_modes,n_msh_elts)
   complex(8) v_beta(n_modes)
   complex(8), dimension(n_modes) :: m_energy
   double precision k_0
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg


   type(NBError) :: nberr


   !f2py intent(in) k_0, n_modes, n_msh_elts, n_msh_pts
   !f2py intent(in) P2_NODES_PER_EL, m_elnd_to_mshpt
   !f2py intent(in) x, v_beta, soln_em_e
   !
   !f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
   !f2py depend(x) n_msh_pts
   !f2py depend(v_beta) n_modes
   !f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_elts
   !
   !f2py intent(out) m_energy


   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call EM_mode_energy_int_impl (k_0, n_modes, n_msh_elts, n_msh_pts,&
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_energy, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine EM_mode_energy_int


subroutine em_mode_act_energy_quadrature (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material, &
   v_refindex, soln_em_e, m_energy, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators

   integer(8) n_modes, n_msh_elts, n_msh_pts
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   integer(8) n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   complex(8) v_refindex(n_elt_mats)
   complex(8) soln_em_e(3,N_DOF_PER_EL,n_modes,n_msh_elts)
   complex(8), dimension(n_modes), intent(out) :: m_energy
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


   type(NBError) nberr

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(v_refindex) n_elt_mats
!f2py depend(v_elt_material) n_msh_elts

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call em_mode_act_energy_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, n_elt_mats, v_elt_material, &
   v_refindex, soln_em_e, m_energy, nberr)

   call finalise_pyerror(nberr, errco, emsg)


end subroutine em_mode_act_energy_quadrature


subroutine em_mode_power_sz_analytic (k_0, n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_power, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators

   double precision k_0      !  k_0 = 2 pi / lambda, where lambda in meters.

   integer(8) n_modes, n_msh_elts, n_msh_pts
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_em_e(3,N_DOF_PER_EL,n_modes,n_msh_elts)
   complex(8) v_beta(n_modes)
   complex(8), dimension(n_modes) :: m_power

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   type(NBError) nberr

!
!f2py intent(in) k_0, n_modes, n_msh_elts, n_msh_pts
!f2py intent(in) P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) x, v_beta, soln_em_e
!
!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(x) n_msh_pts
!f2py depend(v_beta) n_modes
!f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_elts
!
!f2py intent(out) m_power

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call em_mode_power_sz_analytic_impl (k_0, n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_power, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine em_mode_power_sz_analytic



subroutine h_mode_field_ez (k_0, n_modes, n_msh_elts, n_msh_pts, m_elnd_to_mshpt, &
 v_mshpt_xy, v_beta, soln_k1, soln_H1)

   use numbatmod
   integer(8) n_modes, n_msh_elts, n_msh_pts
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_k1(3,P2_NODES_PER_EL+7,n_modes,n_msh_elts)
   complex(8) soln_H1(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8) v_beta(n_modes)
   double precision k_0


   !f2py intent(in) k_0, n_modes, n_msh_elts, n_msh_pts
   !f2py intent(in) P2_NODES_PER_EL, m_elnd_to_mshpt
   !f2py intent(in) x, v_beta, soln_k1
   !
   !f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
   !f2py depend(x) n_msh_pts
   !f2py depend(v_beta) n_modes
   !f2py depend(soln_k1) P2_NODES_PER_EL, n_modes, n_msh_elts
   !
   !f2py intent(out) soln_H1





   !call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call h_mode_field_ez_impl (k_0, n_modes, n_msh_elts, n_msh_pts, m_elnd_to_mshpt, &
 v_mshpt_xy, v_beta, soln_k1, soln_H1)

   !call finalise_pyerror(nberr, errco, emsg)


end subroutine h_mode_field_ez



subroutine moving_boundary (nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac, &
   n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, &
   n_elt_mats, v_elt_material, typ_select_in, typ_select_out, &
   soln_em_p, soln_em_s, soln_ac_u, v_eps_rel, Q_MB, errco, emsg)

   use numbatmod
   integer(8) nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats

   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision v_mshpt_xy(2, n_msh_pts)

   integer(8) typ_select_in, typ_select_out

   complex(8) soln_em_p(3, P2_NODES_PER_EL, nval_em_p, n_msh_elts)
   complex(8) soln_em_s(3, P2_NODES_PER_EL, nval_em_s, n_msh_elts)
   complex(8) soln_ac_u(3, P2_NODES_PER_EL, nval_ac, n_msh_elts)

   complex(8) v_eps_rel(n_elt_mats)
   complex(8), intent(out) :: Q_MB(nval_em_p, nval_em_s, nval_ac)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg


   !f2py intent(in) n_msh_elts, n_msh_pts
   !f2py intent(in) nval_em_p, nval_em_s, nval_ac
   !f2py intent(in) ival_p, ival_s, ival_ac, n_elt_mats
   !f2py intent(in) n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt, debug
   !f2py intent(in) v_elt_material, x, soln_em_p, soln_em_s, soln_ac_u
   !f2py intent(in) typ_select_in, typ_select_out, v_eps_rel, debug

   !f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
   !f2py depend(v_elt_material) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, n_msh_elts
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, n_msh_elts
   !f2py depend(soln_ac_u) P2_NODES_PER_EL, nval_ac, n_msh_elts
   !f2py depend(v_eps_rel) n_elt_mats
   !
   !fo2py intent(out) Q_MB

   type(NBError) nberr

   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call  moving_boundary_impl (nval_em_p, nval_em_s, nval_ac, ival_p, ival_s, ival_ac, &
   n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, &
   n_elt_mats, v_elt_material, typ_select_in, typ_select_out, &
   soln_em_p, soln_em_s, soln_ac_u, v_eps_rel, Q_MB, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine moving_boundary

subroutine photoelastic_int_common (is_curvilinear, nval_em_p, nval_em_s, nval_ac_u, ival_p, &
   ival_s, ival_ac, &
   n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, &
   n_elt_mats, v_elt_material, p_tensor, beta_ac, soln_em_p, soln_em_s, soln_ac_u,&
   v_eps_rel, Q_PE, errco, emsg)

   use numbatmod

   integer(8) is_curvilinear
   integer(8) nval_em_p, nval_em_s, nval_ac_u, ival_p, ival_s, ival_ac
   integer(8) n_msh_elts, n_msh_pts, n_elt_mats
   integer(8) v_elt_material(n_msh_elts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   double precision v_mshpt_xy(2,n_msh_pts)
   complex(8) soln_em_p(3,P2_NODES_PER_EL,nval_em_p,n_msh_elts)
   complex(8) soln_em_s(3,P2_NODES_PER_EL,nval_em_s,n_msh_elts)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,nval_ac_u,n_msh_elts)
   complex(8) p_tensor(3,3,3,3,n_elt_mats)

   complex(8) beta_ac
   complex(8) v_eps_rel(n_elt_mats)

   complex(8), intent(out) :: Q_PE(nval_em_p, nval_em_s, nval_ac_u)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   type(NBError) nberr

   !fo2py intent(in) nval_em_p, nval_em_s, nval_ac_u
   !fo2py intent(in) ival_p, ival_s, ival_ac, n_elt_mats
   !fo2py intent(in) n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt, p_tensor, beta_ac, debug
   !fo2py intent(in) v_elt_material, x, soln_em_p, soln_em_s, soln_ac_u, v_eps_rel
   !
   ! Need these dependencies to get f2py calling to work
   !f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
   !f2py depend(v_elt_material) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, n_msh_elts
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, n_msh_elts
   !f2py depend(soln_ac_u) P2_NODES_PER_EL, nval_ac_u, n_msh_elts
   !f2py depend(p_tensor) n_elt_mats
   !f2py depend(v_eps_rel) n_elt_mats
   !
   !fo2py intent(out) Q_PE
   !fo2py intent(out) errco
   !fo2py intent(out) emsg






   call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call photoelastic_int_common_impl (is_curvilinear, nval_em_p, nval_em_s, nval_ac_u, ival_p, &
   ival_s, ival_ac, &
   n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, &
   n_elt_mats, v_elt_material, p_tensor, beta_ac, soln_em_p, soln_em_s, soln_ac_u,&
   v_eps_rel, Q_PE, nberr)

   call finalise_pyerror(nberr, errco, emsg)

end subroutine photoelastic_int_common



subroutine array_size (n_msh_pts, n_msh_elts, n_modes, &
   int_size, cmplx_size, real_size, n_ddl, errco, emsg)

   use numbatmod

      integer(8), intent(in) :: n_msh_pts, n_msh_elts, n_modes

   integer(8), intent(out) :: int_size, cmplx_size, real_size, n_ddl

   integer(8), intent(out):: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

     type(NBError) nberr

     call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call array_size_impl(n_msh_pts, n_msh_elts, n_modes, &
   int_size, cmplx_size, real_size, n_ddl, nberr)

   call finalise_pyerror(nberr, errco, emsg)



end subroutine



subroutine em_mode_power_sz_quadrature (k_0, n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_power, errco, emsg)

   use numbatmod

   double precision k_0  ! k_0 = 2 pi / lambda, where lambda in meters.

   integer(8) n_modes, n_msh_elts, n_msh_pts
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision v_mshpt_xy(2, n_msh_pts)
   complex(8) soln_em_e(3, N_DOF_PER_EL, n_modes, n_msh_elts)
   complex(8) beta
   complex(8) v_beta(n_modes)
   complex(8), dimension(n_modes) :: m_power

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   type(NBError) nberr

   !f2py intent(in) k_0, n_modes, n_msh_elts, n_msh_pts
!f2py intent(in) P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) x, v_beta, soln_em_e
!
!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(x) n_msh_pts
!f2py depend(v_beta) n_modes
!f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_elts
!
!f2py intent(out) m_power


     call initialise_pyerror(nberr) ! use these to keep f2py quiet. It doesn't like method calls

   call em_mode_power_sz_quadrature_impl (k_0, n_modes, n_msh_elts, n_msh_pts, &
   m_elnd_to_mshpt, v_mshpt_xy, v_beta, soln_em_e, m_power, nberr)

   call finalise_pyerror(nberr, errco, emsg)



end subroutine em_mode_power_sz_quadrature