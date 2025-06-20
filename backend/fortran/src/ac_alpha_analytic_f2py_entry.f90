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
