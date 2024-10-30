#include "numbat_decl.h"

! Calculate the v_power_Sz integral of an AC mode with itself using
! numerical quadrature.

! P_z = Re \int_A \dxdy (-2 i \Omega) c_zjkl u_j^* d_k u_l
!
subroutine ac_mode_power_quadrature (n_modes, n_msh_el, n_msh_pts, &
   v_nd_xy, elnd_to_mshpt,  &
   n_elt_mats, v_el_material, stiffC_zjkl, &
   q_AC, Omega_AC, soln_ac_u, v_power_Sz_r)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, n_msh_el, n_msh_pts
   double precision v_nd_xy(2,n_msh_pts)
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   integer(8) n_elt_mats, v_el_material(n_msh_el)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC
   double precision, dimension(n_modes) :: v_power_Sz_r

   complex(8) stiffC_zjkl(3,3,3,n_elt_mats)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! Locals
   type(NBError) nberr
   complex(8), dimension(n_modes) :: v_power_Sz
   double precision el_nds_xy(2,P2_NODES_PER_EL)

   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar, v_pow, z_tmp1
   integer(8) typ_e, i_el, 
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

!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, stiffC_zjkl, q_AC
!f2py intent(in) soln_ac_u, Omega_AC

!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes

!f2py intent(out) v_power_Sz_r



   errco=0
   call nberr%reset()

   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)


   v_power_Sz = C_ZERO

   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)

      call frontend%nodes_at_el(i_el, el_nds_xy)

      call basfuncs%set_affine_for_elt(el_nds_xy, nberr)
      RET_ON_NBERR(nberr)

      is_curved = frontend%elt_is_curved()
      if (is_curved) n_curved = n_curved + 1


      bas_ovrlp = D_ZERO


      ! For each quadrature point evaluate v_power_Sz of Lagrange polynomials
      ! or derivative of Lagrange polynomials
      do iq=1,quadint%n_quad

         call quadint%get_quad_point(iq, t_xy, t_quadwt)
         RET_ON_NBERR_UNFOLD(nberr)

         call basfuncs%evaluate_at_position(i_el, t_xy, is_curved, el_nds_xy, nberr)
         RET_ON_NBERR(nberr)

         ! Calculate v_power_Sz of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (fi_eld).

         qwt = t_quadwt * abs(basfuncs%det)

         do bf_j=1,P2_NODES_PER_EL
            do xyz_j=1,3
               ind_j = xyz_j + 3*(bf_j-1)

               !  Gradient of transverse components of basis function
               !  sum_{k=x,y} c_zikl u_i d_k u_l
               do xyz_k=1,2

                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)

                        z_tmp1 = basfuncs%phi_P2_ref(bf_j) * basfuncs%gradt_P2_act(xyz_k,bf_l)
                        bas_ovrlp(ind_j,xyz_k,ind_l) = bas_ovrlp(ind_j,xyz_k,ind_l) &
                           + qwt * stiffC_zjkl(xyz_j, xyz_k, xyz_l, typ_e) * z_tmp1
                     enddo
                  enddo

               enddo

               !   Gradient of longitudinal components of basis function,
               !   which is i*beta*phi because fi_eld is assumed to be of
               !   form e^{i*beta*z} phi.

               !   c_zizl u_i d_z u_l = i q c_zizl u_i u_l

               xyz_k=3

               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)

                     z_tmp1 = basfuncs%phi_P2_ref(bf_j) * basfuncs%phi_P2_ref(bf_l) * C_IM_ONE* q_AC

                     bas_ovrlp(ind_j,xyz_k,ind_l) = bas_ovrlp(ind_j,xyz_k,ind_l) &
                        + qwt *  stiffC_zjkl(xyz_j,xyz_k,xyz_l,typ_e) * z_tmp1
                  enddo
               enddo

            enddo
         enddo
      enddo


      ! Having calculated v_power_Sz of basis functions on element
      ! now multiply by specific fi_eld values for modes of interest.
      do md_i=1,n_modes
         v_pow = D_ZERO

         do bf_j=1,P2_NODES_PER_EL
            do xyz_j=1,3
               ind_j = xyz_j + 3*(bf_j-1)

               Ustar = conjg(soln_ac_u(xyz_j,bf_j,md_i,i_el))

               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)
                     U = soln_ac_u(xyz_l,bf_l,md_i,i_el)

                     do xyz_k=1,3
                        v_pow = v_pow  + Ustar * U *  bas_ovrlp(ind_j,xyz_k,ind_l)
                     enddo

                  enddo

               enddo
            enddo
         enddo
         v_power_Sz(md_i) = v_power_Sz(md_i) + v_pow
      enddo
   enddo

   v_power_Sz_r = real(- 2.0d0 * C_IM_ONE* Omega_AC * v_power_Sz)

end subroutine ac_mode_power_quadrature
