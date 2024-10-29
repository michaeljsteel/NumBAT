#include "numbat_decl.h"

! Calculate the v_power_Sz integral of an AC mode with itself using
! numerical quadrature.

! P_z = Re \int_A \dxdy (-2 i \Omega) c_zjkl u_j^* d_k u_l
!
subroutine AC_mode_power_quadrature (n_modes, n_msh_el, n_msh_pts, &
   elnd_to_mesh, v_el_material, v_nd_xy, &
   n_elt_mats, stiffC_zjkl, beta_AC, Omega_AC, soln_ac_u, &
   debug, v_power_Sz)

   use numbatmod
   integer(8) n_modes, md_i
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el), debug
   integer(8) elnd_to_mesh(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8) beta_AC
   complex(8), dimension(n_modes) :: v_power_Sz
   complex(8) stiffC_zjkl(3,3,3,n_elt_mats)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! Locals

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   complex(8) v_pow
   integer(8) i, j, j1, typ_e
   integer(8) iel, ind_i, xyz_i, xyz_k
   integer(8) bf_j, ind_l, xyz_l
   integer(8) bf_i, ui
   complex(8) z_tmp1
   double precision mat_B(2,2), mat_T(2,2)
!
!   NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   integer(8)  n_curved
   logical is_curved
   complex(8) coeff_1, coeff_2
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)
!
!
!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mesh
!f2py intent(in) v_el_material, x, n_elt_mats, stiffC_zjkl, beta_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC
!
!f2py depend(elnd_to_mesh) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(stiffC_zjkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!
!f2py intent(out) v_power_Sz
!


   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "AC_mode_power_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   do i=1,n_modes
      v_power_Sz(i) = 0.0d0
   enddo

   do iel=1,n_msh_el
      typ_e = v_el_material(iel)

      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mesh(j,iel)
         nod_el_p(j) = j1
         xel(:,j) = v_nd_xy(:,j1)
      enddo

      is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif

      bas_ovrlp = 0.0d0


      ! For each quadrature point evaluate v_power_Sz of Lagrange polynomials
      ! or derivative of Lagrange polynomials
      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)

         !  xx   = coordinate on the reference triangle
         !  xx_g = coordinate on the actual triangle
         !  phi2_list = values of Lagrange polynomials (1-6) at each local node.
         !  grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)

         if (.not. is_curved ) then
            !   Rectilinear element
            call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL,&
            &xx_g, det, mat_B, mat_T, errco, emsg)
         else
            !   Isoparametric element, 2024-06-13 fixed version
            call jacobian_p2_2d(xel, P2_NODES_PER_EL, phi2_list,&
            &grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif

         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "AC_m_en_int: det = 0 : iel, det = ", iel, det
            write(*,*) "AC_m_en_int: Aborting..."
            stop
         endif

!  grad_i  = gradient on the actual triangle
!  grad_i  = Transpose(mat_T)*grad_i0
!  Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
         &grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         coeff_1 = ww * abs(det)

         ! Calculate v_power_Sz of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (field).
         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)

               !  Gradient of transverse components of basis function
               !  sum_{k=x,y} c_zikl u_i d_k u_l
               do xyz_k=1,2

                  do bf_j=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_j-1)

                        z_tmp1 = phi2_list(bf_i) * grad2_mat(xyz_k,bf_j)
                        coeff_2 = stiffC_zjkl(xyz_i,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(ind_i,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_k,ind_l) &
                        + coeff_1 * coeff_2 * z_tmp1
                     enddo
                  enddo

               enddo

               !   Gradient of longitudinal components of basis function,
               !   which is i*beta*phi because field is assumed to be of
               !   form e^{i*beta*z} phi.

               !   c_zizl u_i d_z u_l = i q c_zizl u_i u_l

               xyz_k=3

               do bf_j=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_j-1)

                     z_tmp1 = phi2_list(bf_i)&
                     &* phi2_list(bf_j) * C_IM_ONE* beta_AC
                     coeff_2 = stiffC_zjkl(xyz_i,xyz_k,xyz_l,typ_e)
                     bas_ovrlp(ind_i,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_k,ind_l) &
                     + coeff_1 * coeff_2 * z_tmp1
                  enddo
               enddo

            enddo
         enddo
      enddo

      ! Having calculated v_power_Sz of basis functions on element
      ! now multiply by specific field values for modes of interest.
      do md_i=1,n_modes
         v_pow = D_ZERO

         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)

               Ustar = conjg(soln_ac_u(xyz_i,bf_i,md_i,iel))

               do bf_j=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_j-1)
                     U = soln_ac_u(xyz_l,bf_j,md_i,iel)

                     do xyz_k=1,3
                        v_pow = v_pow  + Ustar * U *  bas_ovrlp(ind_i,xyz_k,ind_l)
                     enddo

                  enddo

               enddo
            enddo
         enddo
         v_power_Sz(md_i) = v_pow
      enddo
   enddo

   v_power_Sz = -2.0 * C_IM_ONE* Omega_AC * v_power_Sz

end subroutine AC_mode_power_quadrature
