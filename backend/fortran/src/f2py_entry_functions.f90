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

   integer(8) n_modes, md_i
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

   integer(8) n_modes, ival
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

   integer(8) n_modes, md_i
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

   integer(8) n_modes, md_i
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



   ! Locals
   type(NBError) nberr
   complex(8), dimension(n_modes) :: v_power_Sz
   double precision el_nds_xy(2,P2_NODES_PER_EL)

   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar, v_pow, z_tmp1
   integer(8) typ_e, i_el
   integer(8) bf_j, ind_j, xyz_j, xyz_k
   integer(8) bf_l, ind_l, xyz_l

   integer(8)  iq, md_i
   double precision t_xy(2)
   integer(8)  n_curved
   logical is_curved
   double precision qwt

   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs

   double precision t_quadwt

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
   complex(8) beta1
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
   complex(8) beta, t_power
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