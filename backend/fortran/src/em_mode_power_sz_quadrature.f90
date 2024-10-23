#include "numbat_decl.h"

! Calculate the EM mode power using numerical quadrature for the basis functions.
!
! Sturmberg: Eq. (6)
!
!   S_z = 2 Re[\zhat \dot \int dx dy E^* \cross H] =  2 Re[ zhat \dot  \int dx dy E_t^* \cross H_t]
!       = 2/(\mu_0\omega) Re[E_t^* \dot (\beta E_t + i \nabla_t E_z) ]

subroutine em_mode_power_sz_quadrature (k_0, n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mesh, v_nd_xy, v_beta, soln_em_e, m_power, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators

   double precision k_0  ! k_0 = 2 pi / lambda, where lambda in meters.

   integer(8) n_modes, n_msh_el, n_msh_pts
   integer(8) elnd_to_mesh(P2_NODES_PER_EL, n_msh_el)
   double precision v_nd_xy(2, n_msh_pts)
   complex(8) soln_em_e(3, P2_NODES_PER_EL+7, n_modes, n_msh_el)
   complex(8) beta
   complex(8) v_beta(n_modes)
   complex(8), dimension(n_modes) :: m_power

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) ::  emsg

   ! Local variables
   type(NBError) nberr

   ! Coeffs of transverse components on P2 functions, and longitudinal component on P3 functions
   !   [phi1_x, phi1_y, phi2_x, phi2_y, .., phi6_x, phi6_y, psi1, psi2, .., psi10]
   complex(8) vec_r(2*P2_NODES_PER_EL+10)

   ! Coeffs of transverse components on P2 functions
   !   [phi1_x, phi1_y, phi2_x, phi2_y, .., phi6_x, phi6_y]
   complex(8) vec_l(2*P2_NODES_PER_EL)
   complex(8) vec_1(2*P2_NODES_PER_EL)

   complex(8) bas_ovrlp(2*P2_NODES_PER_EL, 2*P2_NODES_PER_EL+10)

   integer(8) i, j, iq
   integer(8) i_el, ival
   integer(8) nd_j, ind_j, xy_j
   integer(8) nd_i, ind_i, xy_i
   integer(8)  n_curved, offset
   logical is_curved
   double precision nds_xy(2,P2_NODES_PER_EL)
   double precision vec_phi_j(2), vec_phi_i(2)
   complex(8) z_tmp1

   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend

   logical do_P3

   double precision t_quadwt

!f2py intent(in) k_0, n_modes, n_msh_el, n_msh_pts
!f2py intent(in) P2_NODES_PER_EL, elnd_to_mesh
!f2py intent(in) x, v_beta, soln_em_e
!
!f2py depend(elnd_to_mesh) P2_NODES_PER_EL, n_msh_el
!f2py depend(x) n_msh_pts
!f2py depend(v_beta) n_modes
!f2py depend(soln_em_e) P2_NODES_PER_EL, n_modes, n_msh_el
!
!f2py intent(out) m_power
!

   errco = 0
   emsg = ""
   call nberr%reset()

   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mesh, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)

   m_power = D_ZERO

   n_curved = 0
   do_P3 = .true.

   do i_el=1,n_msh_el

      call frontend%nodes_at_el(i_el, nds_xy)

      is_curved = frontend%elt_is_curved()
      if (is_curved) n_curved = n_curved + 1


      bas_ovrlp = D_ZERO

      do iq=1,quadint%n_quad

         call quadint%build_transforms_at(iq, nds_xy, is_curved, do_P3, nberr)
         RET_ON_NBERR_UNFOLD(nberr)

         t_quadwt = quadint%get_current_quadweight()


         do nd_i=1,P2_NODES_PER_EL
            do xy_i=1,2
               ind_i = xy_i + 2*(nd_i-1)

               ! Get basis vector nd_i for the polarisation xy_i
               vec_phi_i = D_ZERO
               vec_phi_i(xy_i) = quadint%phi_P2_ref(nd_i)

               do nd_j=1,P2_NODES_PER_EL
                  do xy_j=1,2
                     ind_j = xy_j + 2*(nd_j-1)

                     ! Get basis vector nd_j for the polarisation xy_j
                     vec_phi_j = D_ZERO
                     vec_phi_j(xy_j) = quadint%phi_P2_ref(nd_j)

                     ! overlaps are diagonal in xx xy \\ yx yy  quadrants
                     z_tmp1 = t_quadwt*(vec_phi_i(1)*vec_phi_j(1) + vec_phi_i(2)*vec_phi_j(2))
                     bas_ovrlp(ind_i,ind_j) = bas_ovrlp(ind_i,ind_j) + z_tmp1
                  enddo
               enddo

               offset = 2*P2_NODES_PER_EL
               do nd_j=1,P3_NODES_PER_EL
                  xy_j = 3
                  ind_j = nd_j + offset

                  ! v_phi_j =  Grad_t psi_nd_j
                  vec_phi_j = quadint%gradt_P3_act(:, nd_j)

                  z_tmp1 = t_quadwt * (vec_phi_i(1)*vec_phi_j(1) + vec_phi_i(2)*vec_phi_j(2))
                  bas_ovrlp(ind_i,ind_j) = bas_ovrlp(ind_i,ind_j) + z_tmp1
               enddo

            enddo
         enddo
      enddo

      bas_ovrlp = bas_ovrlp


      do ival=1,n_modes

         beta = v_beta(ival)

         ! v_l = e_t^*
         do nd_i=1,P2_NODES_PER_EL
            do xy_j=1,2
               ind_i = xy_j + 2*(nd_i-1)
               vec_l(ind_i) = conjg(soln_em_e(xy_j, nd_i, ival, i_el))
            enddo
         enddo

         ! v_r = beta e_t ...
         do nd_i=1,P2_NODES_PER_EL
            do xy_j=1,2
               ind_j = xy_j + 2*(nd_i-1)
               vec_r(ind_j) = soln_em_e(xy_j, nd_i, ival, i_el) * beta
            enddo
         enddo

         ! v_r += i nabla e_z
         offset = 2*P2_NODES_PER_EL
         ! The longitudinal component at the vertices (P3 elements) = -i Ez
         do nd_i=P3_VERT_1, P3_VERT_3
            ind_j = nd_i + offset
            vec_r(ind_j) =  soln_em_e(3, nd_i, ival, i_el) * C_IM_ONE
         enddo

         do nd_i=P3_EDGE_LO, P3_INTERIOR
            ind_j = offset + nd_i

            vec_r(ind_j) =  soln_em_e(3, nd_i+P3_VERT_3, ival, i_el) * C_IM_ONE
         enddo

         ! Calculate  vec_l^t M_ij vec_r

         ! MH = M_ij vec_r
         do i=1,2*P2_NODES_PER_EL
            vec_1(i) = 0.0d0
            do j=1,2*P2_NODES_PER_EL + P3_NODES_PER_EL
               vec_1(i) = vec_1(i) +  bas_ovrlp(i,j) * vec_r(j)
            enddo
         enddo

         ! S_z = vec_r^t MH
         do i=1,2*P2_NODES_PER_EL
            m_power(ival) = m_power(ival) + vec_l(i) * vec_1(i)
         enddo
      enddo
   enddo

   m_power = m_power/ (k_0 * SI_C_SPEED * SI_MU_0)

end subroutine em_mode_power_sz_quadrature
