#include "numbat_decl.h"

!  Calculate the v_alpha integral of an AC mode with itself using
!  numerical quadrature.

! \alpha = \Omega^2/Energy_aC \int  eta_ijkl d_i u_j^* d_k u_l


subroutine ac_alpha_quadrature_impl (n_modes, n_msh_elts, n_msh_pts, &
   v_mshpt_xy, m_elnd_to_mshpt, n_elt_mats, v_elt_material, eta_ijkl, &
   q_AC, Omega_AC, soln_ac_u, &
   v_ac_mode_energy, v_alpha_r, nberr)


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
   type(NBError) nberr

   ! Locals
   complex(8), dimension(n_modes) :: v_alpha
   double precision el_nds_xy(2,P2_NODES_PER_EL)

   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) typ_e
   integer(8) i_el, ind_i, xyz_i, xyz_k
   integer(8) bf_l, ind_l, xyz_l, xyz_j
   integer(8) bf_i, iq
   complex(8) z_tmp1
   double precision t_xy(2)

!     NQUAD: The number of quadrature points used in each element.
   integer(8)  n_curved
   logical is_curved
   complex(8) t_eta

   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs

   double precision t_quadwt, t_qwt


   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, nberr)
   RET_ON_NBERR(nberr)

   v_alpha = C_ZERO

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)

      call frontend%nodes_at_el(i_el, el_nds_xy)

      call basfuncs%set_affine_for_elt(el_nds_xy, nberr)
      RET_ON_NBERR(nberr)

      is_curved = frontend%elt_is_curved()

      if (is_curved) n_curved = n_curved + 1
      bas_ovrlp  = D_ZERO

      ! For each quadrature point evaluate v_alpha of Lagrange polynomials
      ! or derivative of Lagrange polynomials
      do iq=1,quadint%n_quad

         call quadint%get_quad_point(iq, t_xy, t_quadwt)
         RET_ON_NBERR(nberr)

         call basfuncs%evaluate_at_position(i_el, t_xy, is_curved, el_nds_xy, nberr)
         RET_ON_NBERR(nberr)

         t_qwt = t_quadwt * abs(basfuncs%det)


         ! Calculate v_alpha of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function ().
         ! TODO: swap meaning of i and j to fit the equation at top.
         !        current version is correct but confusing. ac_alpha_analytic has made the switch
         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)

               !         Gradient of transverse components of basis function
               do xyz_j=1,2
                  do xyz_k=1,2
                     do bf_l=1,P2_NODES_PER_EL
                        do xyz_l=1,3
                           ind_l = xyz_l + 3*(bf_l-1)

                           z_tmp1 = basfuncs%gradt_P2_act(xyz_j,bf_i) * basfuncs%gradt_P2_act(xyz_k,bf_l)
                           t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)

                           bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                              + t_qwt * t_eta * z_tmp1
                        enddo
                     enddo
                  enddo

                  !  Gradient of longitudinal components of basis function,
                  !  which is i*beta*phi because  is assumed to be of form e^{i*beta*z} phi.
                  xyz_k=3
                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)

                        z_tmp1 = basfuncs%gradt_P2_act(xyz_j,bf_i) * basfuncs%phi_P2_ref(bf_l) * C_IM_ONE* q_AC
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                           + t_qwt * t_eta * z_tmp1
                     enddo
                  enddo
               enddo

               xyz_j=3
               do xyz_k=1,2
                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)
                        z_tmp1 = basfuncs%phi_P2_ref(bf_i) * (-C_IM_ONE* q_AC) * basfuncs%gradt_P2_act(xyz_k,bf_l)
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) =  bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                           + t_qwt * t_eta * z_tmp1
                     enddo
                  enddo
               enddo

               xyz_k=3
               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)
                     z_tmp1 = q_AC**2 * basfuncs%phi_P2_ref(bf_i) * basfuncs%phi_P2_ref(bf_l)
                     t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                     bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                        + t_qwt * t_eta * z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! Having calculated v_alpha of basis functions on element
      ! now multiply by specific  values for modes of interest.
      do ival=1,n_modes

         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)

               Ustar = conjg(soln_ac_u(xyz_i,bf_i,ival,i_el))

               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)

                     U = soln_ac_u(xyz_l,bf_l,ival,i_el)

                     do xyz_j=1,3
                        do xyz_k=1,3
                           z_tmp1 = Ustar * U * bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l)
                           v_alpha(ival) = v_alpha(ival) + z_tmp1
                        enddo
                     enddo

                  enddo
               enddo

            enddo
         enddo

      enddo

   enddo

!    ! Multiply through prefactor
!    do i=1,n_modes
! !         z_tmp1 = -1.0 * Omega_AC(i)**2 / v_ac_mode_energy(i)
! !       Flipped sign as assuming did not do integration by parts - going off CW advice.
!       z_tmp1 = Omega_AC(i)**2 / v_ac_mode_energy(i)
!       v_alpha(i) = z_tmp1 * v_alpha(i)
!    enddo

   v_alpha_r = real(v_alpha * Omega_AC**2 / v_ac_mode_energy)

end subroutine ac_alpha_quadrature_impl
