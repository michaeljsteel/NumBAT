! Calculate the overlap integral of an AC mode with itself using
! numerical quadrature.
!
subroutine AC_mode_power_int (nval,&
&nel, npt, nnodes, elnd_to_mesh, type_el, x,&
&nb_typ_el, c_tensor_z, beta_AC, Omega_AC, soln_AC,&
&debug, overlap)
!
   use numbatmod
   integer(8) nval, ival
   integer(8) nel, npt, nnodes, nb_typ_el
   integer(8) type_el(nel), debug
   integer(8) elnd_to_mesh(nnodes,nel)
   double precision x(2,npt)
!      complex(8) x(2,npt)
   complex(8) soln_AC(3,nnodes,nval,nel)
   complex(8) Omega_AC(nval)
   complex(8) beta_AC
   complex(8), dimension(nval) :: overlap
   complex(8) c_tensor_z(3,3,3,nb_typ_el)
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

!     Local variables

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) basis_overlap(3*P2_NODES_PER_EL,3,3*P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) i, j, k, l, j1, typ_e
   integer(8) iel, ind_ip, i_eq, k_eq
   integer(8) ltest, ind_lp, l_eq
   integer(8) itrial, ui
   complex(8) z_tmp1
   double precision mat_B(2,2), mat_T(2,2)
!
!     NQUAD: The number of quadrature points used in each element.
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
!f2py intent(in) nval, nel, npt, nnodes, elnd_to_mesh
!f2py intent(in) type_el, x, nb_typ_el, c_tensor_z, beta_AC
!f2py intent(in) soln_AC, debug, Omega_AC
!
!f2py depend(elnd_to_mesh) nnodes, nel
!f2py depend(type_el) npt
!f2py depend(x) npt
!f2py depend(soln_AC) nnodes, nval, nel
!f2py depend(c_tensor_z) nb_typ_el
!f2py depend(Omega_AC) nval
!
!f2py intent(out) overlap
!
!
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
!
   ui = stdout
!
   if ( nnodes .ne. 6 ) then
      write(ui,*) "AC_mode_power_int: problem nnodes = ", nnodes
      write(ui,*) "AC_mode_power_int: nnodes should be equal to 6 !"
      write(ui,*) "AC_mode_power_int: Aborting..."
      stop
   endif
!
   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "AC_mode_power_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif
!
   do i=1,nval
      overlap(i) = 0.0d0
   enddo
!
!ccccccccccc
! Loop over elements - start
!ccccccccccc
   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,nnodes
         j1 = elnd_to_mesh(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      is_curved = log_is_curved_elem_tri (nnodes, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif
!ccccccccc
      do i=1,3*nnodes
         do k=1,3
            do l=1,3*nnodes
               basis_overlap(i,k,l) = 0.0d0
            enddo
         enddo
      enddo
!ccccccccc
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
         if (.not. is_curved ) then
!           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,&
            &xx_g, det, mat_B, mat_T, errco, emsg)
         else
!           Isoparametric element, 2024-06-13 fixed version
            call jacobian_p2_2d(xel, nnodes, phi2_list,&
            &grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif
         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "AC_m_en_int: det = 0 : iel, det = ", iel, det
            write(*,*) "AC_m_en_int: Aborting..."
            stop
         endif
!          grad_i  = gradient on the actual triangle
!          grad_i  = Transpose(mat_T)*grad_i0
!          Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
         &grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         coeff_1 = ww * abs(det)
! Calculate overlap of basis functions at quadrature point,
! which is a superposition of P2 polynomials for each function (field).
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
!             Gradient of transverse components of basis function
               do k_eq=1,2
                  do ltest=1,P2_NODES_PER_EL
                     do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
                        z_tmp1 = phi2_list(itrial) * grad2_mat(k_eq,ltest)
                        coeff_2 = c_tensor_z(i_eq,k_eq,l_eq,typ_e)
                        basis_overlap(ind_ip,k_eq,ind_lp) =&
                        &basis_overlap(ind_ip,k_eq,ind_lp)&
                        &+ coeff_1 * coeff_2 * z_tmp1
                     enddo
                  enddo
               enddo
!             Gradient of longitudinal components of basis function,
!             which is i*beta*phi because field is assumed to be of
!             form e^{i*beta*z} phi.
               k_eq=3
               do ltest=1,P2_NODES_PER_EL
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     z_tmp1 = phi2_list(itrial)&
                     &* phi2_list(ltest) * C_IM_ONE* beta_AC
                     coeff_2 = c_tensor_z(i_eq,k_eq,l_eq,typ_e)
                     basis_overlap(ind_ip,k_eq,ind_lp) =&
                     &basis_overlap(ind_ip,k_eq,ind_lp)&
                     &+ coeff_1 * coeff_2 * z_tmp1
                  enddo
               enddo
            enddo
         enddo
      enddo
!ccccccccc
! Having calculated overlap of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,nval
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               ind_ip = i_eq + 3*(itrial-1)
               Ustar = conjg(soln_AC(i_eq,itrial,ival,iel))
               do ltest=1,P2_NODES_PER_EL
                  do l_eq=1,3
                     ind_lp = l_eq + 3*(ltest-1)
                     U = soln_AC(l_eq,ltest,ival,iel)
                     do k_eq=1,3
                        z_tmp1 = basis_overlap(ind_ip,k_eq,ind_lp)
                        overlap(ival) = overlap(ival) + Ustar * U * z_tmp1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!ccccccccccc
! Loop over elements - end
!ccccccccccc
   enddo
! Multiply through prefactor
   do i=1,nval
      overlap(i) = -2.0 * C_IM_ONE* Omega_AC(i) * overlap(i)
   enddo
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end subroutine AC_mode_power_int
