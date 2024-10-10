!  Calculate the overlap integral of an AC mode with itself using
!  numerical quadrature.

#include "numbat_decl.h"

subroutine AC_alpha_int (nval,&
nel, npt, nnodes, table_nod, type_el, x,&
nb_typ_el, eta_tensor, beta_AC, Omega_AC, soln_AC,&
AC_mode_energy_elastic, debug, overlap, errco, emsg)


   use numbatmod

   integer(8) nval, ival
   integer(8) nel, npt, nnodes, nb_typ_el
   integer(8) type_el(nel), debug
   integer(8) table_nod(nnodes,nel)
   double precision x(2,npt)
!       complex(8) x(2,npt)
   complex(8) soln_AC(3,nnodes,nval,nel)
   complex(8) Omega_AC(nval)
   complex(8) beta_AC, AC_mode_energy_elastic(nval)
   complex(8), dimension(nval) :: overlap
   complex(8) eta_tensor(3,3,3,3,nb_typ_el)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

!     Local variables
   integer(8) nnodes0
   parameter (nnodes0 = 6)
   integer(8) nod_el_p(nnodes0)
   double precision xel(2,nnodes0)
   complex(8) basis_overlap(3*nnodes0,3,3,3*nnodes0)
   complex(8) U, Ustar
   integer(8) i, j, k, l, j1, typ_e
   integer(8) iel, ind_ip, i_eq, k_eq
   integer(8) ltest, ind_lp, l_eq, j_eq
   integer(8) itrial, ui
   complex(8) z_tmp1
   double precision mat_B(2,2), mat_T(2,2)

!     NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   integer(8)  n_curved
   logical info_curved
   complex(8) coeff_1, coeff_2
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)
!
!
!f2py intent(in) nval, nel, npt, nnodes, table_nod
!f2py intent(in) type_el, x, nb_typ_el, eta_tensor, beta_AC
!f2py intent(in) soln_AC, debug, Omega_AC, AC_mode_energy_elastic
!
!f2py depend(table_nod) nnodes, nel
!f2py depend(type_el) npt
!f2py depend(x) npt
!f2py depend(soln_AC) nnodes, nval, nel
!f2py depend(eta_tensor) nb_typ_el
!f2py depend(Omega_AC) nval
!f2py depend(AC_mode_energy_elastic) nval
!
!f2py intent(out) overlap
!
!

   ui = stdout

!
   if ( nnodes .ne. 6 ) then
      write(ui,*) "AC_alpha_int: problem nnodes = ", nnodes
      write(ui,*) "AC_alpha_int: nnodes should be equal to 6 !"
      write(ui,*) "AC_alpha_int: Aborting..."
      stop
   endif
!
   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "AC_alpha_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif
!
   do i=1,nval
      overlap(i) = 0.0d0
   enddo


   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = table_nod(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      info_curved = log_is_curved_elem_tri (nnodes, xel)
      if (info_curved) then
         n_curved = n_curved + 1
      endif


      do i=1,3*nnodes
         do j=1,3
            do k=1,3
               do l=1,3*nnodes
                  basis_overlap(i,j,k,l) = 0.0d0
               enddo
            enddo
         enddo
      enddo

! For each quadrature point evaluate overlap of Lagrange polynomials
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
         if (.not. info_curved) then
!           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,&
            xx_g, det, mat_B, mat_T, errco, emsg)
            RETONERROR(errco)
         else
!           Isoparametric element!  2024-06-13 fixed version
            call jacobian_p2_2d(xel, nnodes, phi2_list,&
            grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
            RETONERROR(errco)
         endif


           grad_i  = gradient on the actual triangle
!          grad_i  = Transpose(mat_T)*grad_i0
!          Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
         grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         coeff_1 = ww * abs(det)
! Calculate overlap of basis functions at quadrature point,
! which is a superposition of P2 polynomials for each function (field).
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
!         Gradient of transverse components of basis function
               do j_eq=1,2
                  do k_eq=1,2
                     do ltest=1,nnodes0
                        do l_eq=1,3
                           ind_lp = l_eq + 3*(ltest-1)
                           z_tmp1 = grad2_mat(j_eq,itrial)&
                           &* grad2_mat(k_eq,ltest)
                           coeff_2 = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                           basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =&
                           &basis_overlap(ind_ip,j_eq,k_eq,ind_lp)&
                           &+ coeff_1 * coeff_2 * z_tmp1
                        enddo
                     enddo
                  enddo
!           Gradient of longitudinal components of basis function,
!           which is i*beta*phi because field is assumed to be of
!           form e^{i*beta*z} phi.
                  k_eq=3
                  do ltest=1,nnodes0
                     do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
                        z_tmp1 = grad2_mat(j_eq,itrial)&
                        &* phi2_list(ltest) * C_IM_ONE* beta_AC
                        coeff_2 = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                        basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =&
                        &basis_overlap(ind_ip,j_eq,k_eq,ind_lp)&
                        &+ coeff_1 * coeff_2 * z_tmp1
                     enddo
                  enddo
               enddo

               j_eq=3
               do k_eq=1,2
                  do ltest=1,nnodes0
                     do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
                        z_tmp1 = phi2_list(itrial) * (-C_IM_ONE* beta_AC)&
                        &* grad2_mat(k_eq,ltest)
                        coeff_2 = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                        basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =&
                        &basis_overlap(ind_ip,j_eq,k_eq,ind_lp)&
                        &+ coeff_1 * coeff_2 * z_tmp1
                     enddo
                  enddo
               enddo

               k_eq=3
               do ltest=1,nnodes0
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     z_tmp1 = beta_AC**2 * phi2_list(itrial)&
                     &* phi2_list(ltest)
                     coeff_2 = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                     basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =&
                     &basis_overlap(ind_ip,j_eq,k_eq,ind_lp)&
                     &+ coeff_1 * coeff_2 * z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo

! Having calculated overlap of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,nval
         do itrial=1,nnodes0
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               Ustar = conjg(soln_AC(i_eq,itrial,ival,iel))
               do ltest=1,nnodes0
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     U = soln_AC(l_eq,ltest,ival,iel)
                     do j_eq=1,3
                        do k_eq=1,3
                           z_tmp1 = basis_overlap(ind_ip,j_eq,k_eq,ind_lp)
                           z_tmp1 = Ustar * U * z_tmp1
                           overlap(ival) = overlap(ival) + z_tmp1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

   enddo

   ! Multiply through prefactor
   do i=1,nval
!         z_tmp1 = -1.0 * Omega_AC(i)**2 / AC_mode_energy_elastic(i)
!       Flipped sign as assuming did not do integration by parts - going off CW advice.
      z_tmp1 = Omega_AC(i)**2 / AC_mode_energy_elastic(i)
      overlap(i) = z_tmp1 * overlap(i)
   enddo

end subroutine AC_alpha_int
