#include "numbat_decl.h"

! Calculate the v_power_Sz integral of an AC mode with itself using
! analytic expressions for basis function v_power_Szs on linear elements.

! P_z = Re (-2 i \Omega) \int_A \dxdy  c_zjkl u_j^* d_k u_l


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


   ! Locals

   complex(8), dimension(n_modes):: v_power_Sz
   double precision nds_xy(2, P2_NODES_PER_EL)
   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3*P2_NODES_PER_EL)
   complex(8) U, Ustar, v_pow
   integer(8) typ_e, md_i
   integer(8) i_el, bf_i, ind_i, xyz_i
   integer(8) bf_j, ind_j, xyz_j

   complex(8) stiff_C_IJ(6,6)

   type(NBError) nberr
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs

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


   errco = 0
   call nberr%reset()


   call frontend%init_from_py(n_msh_elts, n_msh_pts, m_elnd_to_mshpt, v_mshpt_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)

   v_power_Sz = 0.0d0

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)
      stiff_C_IJ = stiffC_IJ_el(:,:,typ_e)

      call frontend%nodes_at_el(i_el, nds_xy)

      call basfuncs%set_affine_for_elt(nds_xy, nberr)
      RET_ON_NBERR(nberr)


      call mat_el_powerflow (q_AC, stiff_C_IJ, basfuncs, bas_ovrlp)


      ! Having calculated v_power_Sz of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes
         v_pow = D_ZERO

         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)
               Ustar = conjg(soln_ac_u(xyz_i,bf_i,md_i,i_el))

               do bf_j=1,P2_NODES_PER_EL
                  do xyz_j=1,3
                     ind_j = xyz_j + 3*(bf_j-1)
                     U = soln_ac_u(xyz_j,bf_j,md_i,i_el)

                     v_pow = v_pow + Ustar * U * bas_ovrlp(ind_i,ind_j)
                  enddo
               enddo

            enddo
         enddo
         v_power_Sz(md_i) = v_power_Sz(md_i) + v_pow
      enddo
   enddo


   v_power_Sz_r = real(- 2.0d0 * C_IM_ONE* Omega_AC * v_power_Sz)


end subroutine ac_mode_power_analytic







subroutine AC_mode_power_analytic_old (n_modes, n_msh_elts, n_msh_pts,  &
   v_mshpt_xy, m_elnd_to_mshpt, &
   n_elt_mats, v_elt_material, stiffC_zjkl, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz)

   use numbatmod

   integer(8) n_modes, n_msh_elts, n_msh_pts
   double precision v_mshpt_xy(2,n_msh_pts)
   integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL,n_msh_elts)
   integer(8) n_elt_mats
   integer(8) v_elt_material(n_msh_elts)

   complex(8) stiffC_zjkl(6,6,n_elt_mats)

   complex(8) q_AC
   complex(8) Omega_AC(n_modes)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_elts)
   complex(8), dimension(n_modes) :: v_power_Sz

   ! Locals

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) j, j1, typ_e
   integer(8) i_el, ind_i, xyz_i
   integer(8) bf_j, ind_j, xyz_j
   integer(8) bf_i, md_i

   complex(8) z_tmp1, stiff_C_IJ(6,6)



!f2py intent(in) n_modes, n_msh_elts, n_msh_pts, P2_NODES_PER_EL, m_elnd_to_mshpt
!f2py intent(in) v_elt_material, x, n_elt_mats, stiffC_zjkl, q_AC
!f2py intent(in) soln_ac_u, Omega_AC

!f2py depend(m_elnd_to_mshpt) P2_NODES_PER_EL, n_msh_elts
!f2py depend(v_elt_material) n_msh_pts
!f2py depend(v_mshpt_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_elts
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes

!f2py intent(out) v_power_Sz




   v_power_Sz = 0.0d0

   do i_el=1,n_msh_elts
      typ_e = v_elt_material(i_el)
      do j=1,P2_NODES_PER_EL
         j1 = m_elnd_to_mshpt(j,i_el)
         nod_el_p(j) = j1
         xel(:,j) = v_mshpt_xy(:,j1)
      enddo

      stiff_C_IJ = stiffC_zjkl(:,:,typ_e)

      call mat_el_powerflow_old (xel, q_AC, stiff_C_IJ, bas_ovrlp)


      ! Having calculated v_power_Sz of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes
         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)
               Ustar = conjg(soln_ac_u(xyz_i,bf_i,md_i,i_el))
               do bf_j=1,P2_NODES_PER_EL
                  do xyz_j=1,3
                     ind_j = xyz_j + 3*(bf_j-1)
                     U = soln_ac_u(xyz_j,bf_j,md_i,i_el)
                     z_tmp1 = bas_ovrlp(ind_i,ind_j)
                     v_power_Sz(md_i) = v_power_Sz(md_i) + Ustar * U * z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo


   v_power_Sz = 2.0d0 * C_IM_ONE* Omega_AC * v_power_Sz

end subroutine AC_mode_power_analytic_old
