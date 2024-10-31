#include "numbat_decl.h"

!  Calculate the v_alpha integral of an AC mode with itself using
!  numerical quadrature.

! \alpha = \Omega^2/Energy_aC \int  eta_ijkl d_i u_j^* d_k u_l



subroutine ac_alpha_quadrature (n_modes, n_msh_el, n_msh_pts, &
   v_nd_xy, elnd_to_mshpt, n_elt_mats, v_el_material, eta_ijkl, &
   q_AC, Omega_AC, soln_ac_u, &
   AC_mode_energy_elastic, debug, v_alpha)


   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   integer(8) n_modes, ival
   integer(8) n_msh_el, n_msh_pts,  n_elt_mats
   integer(8) v_el_material(n_msh_el), debug
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)

   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8) q_AC, AC_mode_energy_elastic(n_modes)
   complex(8), dimension(n_modes) :: v_alpha
   complex(8) eta_ijkl(3,3,3,3,n_elt_mats)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

!     Local variables

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   double precision el_nds_xy(2,P2_NODES_PER_EL)

   complex(8) bas_ovrlp(3*P2_NODES_PER_EL,3,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) i, j, k, l, j1, typ_e
   integer(8) i_el, ind_i, xyz_i, xyz_k
   integer(8) bf_l, ind_l, xyz_l, xyz_j
   integer(8) bf_i, ui
   complex(8) z_tmp1
   double precision mat_B(2,2), mat_T(2,2)

!     NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   integer(8)  n_curved
   logical is_curved
   complex(8) t_eta
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)

type(NBError) nberr
   type(QuadIntegrator) quadint
   type(PyFrontEnd) frontend
   type(BasisFunctions) basfuncs

   double precision t_quadwt, t_qwt

!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, eta_ijkl, q_AC
!f2py intent(in) soln_ac_u, debug, Omega_AC, AC_mode_energy_elastic
!
!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(v_nd_xy) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(eta_ijkl) n_elt_mats
!f2py depend(Omega_AC) n_modes
!f2py depend(AC_mode_energy_elastic) n_modes
!
!f2py intent(out) v_alpha


   ui = stdout

   errco=0
   call nberr%reset()

   call quadint%setup_reference_quadratures()

   call frontend%init_from_py(n_msh_el, n_msh_pts, elnd_to_mshpt, v_nd_xy, nberr)
   RET_ON_NBERR_UNFOLD(nberr)



   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "AC_alpha_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   v_alpha = C_ZERO


   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)

      call frontend%nodes_at_el(i_el, el_nds_xy)

      call basfuncs%set_affine_for_elt(el_nds_xy, nberr)
      RET_ON_NBERR(nberr)

      is_curved = frontend%elt_is_curved()
      if (is_curved) n_curved = n_curved + 1


      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mshpt(j,i_el)
         nod_el_p(j) = j1
         xel(1,j) = v_nd_xy(1,j1)
         xel(2,j) = v_nd_xy(2,j1)
      enddo

      ! info_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, xel)
      ! if (info_curved) then
      !    n_curved = n_curved + 1
      ! endif


      bas_ovrlp  = D_ZERO

! For each quadrature point evaluate v_alpha of Lagrange polynomials
! or derivative of Lagrange polynomials
      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
!         xx   = coordinate on the reference triangle
!         xx_g = coordinate on the actual triangle
!         phi2_list = values of Lagrange polynomials (1-6) at each local node.
!         grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)
!
         if (.not. is_curved) then
!           Rectilinear element
            call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL,&
               xx_g, det, mat_B, mat_T, errco, emsg)
            RETONERROR(errco)
         else
!           Isoparametric element!  2024-06-13 fixed version
            call jacobian_p2_2d(xel, P2_NODES_PER_EL, phi2_list,&
               grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
            RETONERROR(errco)
         endif


         ! grad_i  = gradient on the actual triangle
!          grad_i  = Transpose(mat_T)*grad_i0
!          Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
            grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         t_quadwt = ww * abs(det)

         ! Calculate v_alpha of basis functions at quadrature point,
         ! which is a superposition of P2 polynomials for each function (fi_eld).
         do bf_i=1,P2_NODES_PER_EL
            do xyz_i=1,3
               ind_i = xyz_i + 3*(bf_i-1)

               !         Gradient of transverse components of basis function
               do xyz_j=1,2
                  do xyz_k=1,2
                     do bf_l=1,P2_NODES_PER_EL
                        do xyz_l=1,3
                           ind_l = xyz_l + 3*(bf_l-1)

                           z_tmp1 = grad2_mat(xyz_j,bf_i) * grad2_mat(xyz_k,bf_l)
                           t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)

                           bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                              + t_quadwt * t_eta * z_tmp1
                        enddo
                     enddo
                  enddo

                  !           Gradient of longitudinal components of basis function,
                  !           which is i*beta*phi because fi_eld is assumed to be of
                  !           form e^{i*beta*z} phi.
                  xyz_k=3
                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)

                        z_tmp1 = grad2_mat(xyz_j,bf_i) * phi2_list(bf_l) * C_IM_ONE* q_AC
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                           + t_quadwt * t_eta * z_tmp1
                     enddo
                  enddo
               enddo

               xyz_j=3
               do xyz_k=1,2
                  do bf_l=1,P2_NODES_PER_EL
                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(bf_l-1)
                        z_tmp1 = phi2_list(bf_i) * (-C_IM_ONE* q_AC) * grad2_mat(xyz_k,bf_l)
                        t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                        bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) =  bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                           + t_quadwt * t_eta * z_tmp1
                     enddo
                  enddo
               enddo

               xyz_k=3
               do bf_l=1,P2_NODES_PER_EL
                  do xyz_l=1,3
                     ind_l = xyz_l + 3*(bf_l-1)
                     z_tmp1 = q_AC**2 * phi2_list(bf_i) * phi2_list(bf_l)
                     t_eta = eta_ijkl(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
                     bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) = bas_ovrlp(ind_i,xyz_j,xyz_k,ind_l) &
                        + t_quadwt * t_eta * z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo

! Having calculated v_alpha of basis functions on element
! now multiply by specific fi_eld values for modes of interest.
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
! !         z_tmp1 = -1.0 * Omega_AC(i)**2 / AC_mode_energy_elastic(i)
! !       Flipped sign as assuming did not do integration by parts - going off CW advice.
!       z_tmp1 = Omega_AC(i)**2 / AC_mode_energy_elastic(i)
!       v_alpha(i) = z_tmp1 * v_alpha(i)
!    enddo

   v_alpha = v_alpha * Omega_AC**2 / AC_mode_energy_elastic

end subroutine ac_alpha_quadrature
