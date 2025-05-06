#include "numbat_decl.h"

! Calculate the EM mode power using analytic expressions for the basis functions.
!
! Sturmberg: Eq. (6)
!
!   P_z = 2 Re[\zhat \dot \int dx dy E^* \cross H] =  2 Re[ zhat \dot  \int dx dy E_t^* \cross H_t]
!

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


   ! Locals

   double precision nds_xy(2,P2_NODES_PER_EL)

   complex(8) t_sol_E(3,P2_NODES_PER_EL)
   complex(8) t_sol_H(3,P2_NODES_PER_EL)
   complex(8) t_sol_Ez_P3(P3_NODES_PER_EL)

   !  P3 Ez-field
   double precision m_int_p2_p2(P2_NODES_PER_EL, P2_NODES_PER_EL)
   integer(8) i_el, ival
   integer(8) nd_i, nd_j
   complex(8) vec_Es(3), vec_H(3)
   complex(8) t_Pz

   type(AnalyticIntegrator) integrator   ! TODO: replace with BasicFunctions object
   type(PyFrontEnd) frontend
   integer(8) ilo, ihi, off

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


   call nberr%reset()


   call frontend%init_from_py(n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)

   do ival=1,n_modes
      t_power = D_ZERO
      beta = v_beta(ival)

      do i_el=1,n_msh_elts

         call frontend%nodes_at_el(i_el, nds_xy)

         call integrator%build_transforms_at(nds_xy, nberr)
         RET_ON_NBERR_UNFOLD(nberr)

         !  The matrix m_int_p2_p2 contains the overlap integrals between the P2-polynomial basis functions
         call find_overlaps_p2_p2(m_int_p2_p2, integrator%det)

         ! Need the Et and Ht fields at the P2 nodes
         ! Getting the Ht fields requires the Ez field which requires the P3 longitudinal solutions

         !  The components (E_x,E_y) of the mode ival
         !  The component E_z of the mode ival.
         ! The FEM code uses the scaling: E_z = C_IM_ONE* beta * \hat{E}_z
         t_sol_E = soln_em_e(:, 1:P2_NODES_PER_EL, ival, i_el)

         !  E_z-field: The longitudinal component at the P2 vertices, which are also P3 elements
         ilo = P3_VERT_1
         ihi = P3_VERT_3
         t_sol_Ez_P3(ilo:ihi) = soln_em_e(3, ilo:ihi, ival, i_el)

         !  The longitudinal component at the edge nodes and interior node (P3 elements)
         ilo = P3_EDGE_LO
         ihi = P3_INTERIOR
         off = P3_VERT_3
         t_sol_Ez_P3(ilo:ihi) = soln_em_e(3, ilo+off:ihi+off, ival, i_el)

         !ilo = P2_NODES_PER_EL+1
         !ihi = P2_NODES_PER_EL+P3_NODES_PER_EL-3
         !t_sol_Ez_P3(4:P3_NODES_PER_EL) = soln_em_e(3, ilo:ihi, ival, i_el)

         call get_H_field_p3 (k_0, beta, integrator%mat_T, t_sol_E, t_sol_Ez_P3, t_sol_H)


         do nd_i=1,P2_NODES_PER_EL
            vec_Es = t_sol_E(:, nd_i)

            do nd_j=1,P2_NODES_PER_EL
               vec_H = t_sol_H(:, nd_j)

               !   Cross-product Z.(E^* X H) of E^*=vec_Es and H=vec_H
               !TODO: doesn't seem to be conjugating E field. Doesn't matter since transverse fields are real
               t_Pz = vec_Es(1) * vec_H(2) - vec_Es(2) * vec_H(1)
               t_power = t_power + t_Pz * m_int_p2_p2(nd_i, nd_j)
            enddo
         enddo

      enddo

      m_power(ival) = t_power
   enddo

end subroutine em_mode_power_sz_analytic
