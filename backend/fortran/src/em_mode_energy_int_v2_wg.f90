! Calculate the overlap integral of an EM mode with itself using
! analytic expressions for basis function overlaps on linear elements.
!
subroutine EM_mode_energy_int_v2_wg (k_0, nval, nel, npt,&
&nnodes_P2, m_elnd_to_mshpt,&
&x, betas, soln_k1, type_el, overlap)
!
!     k_0 = 2 pi / lambda, where lambda in meters.
!
   use numbatmod
   integer(8) nval, nel, npt, nnodes_P2
   integer(8) m_elnd_to_mshpt(nnodes_P2,nel)
   double precision x(2,npt)
!      complex(8) x(2,npt)
   complex(8) soln_k1(3,nnodes_P2+7,nval,nel)
   complex(8) beta1, overlap1
   complex(8) betas(nval)
   complex(8), dimension(nval) :: overlap
   double precision k_0

!     Local variables

   integer(8) nnodes_P2_0, nnodes_P3_0
   parameter (nnodes_P2_0 = 6)
   parameter ( nnodes_P3_0 = 10)
   integer(8) nod_el_p(nnodes_P2_0)
   integer(8) type_el(nel), typ_e
   double precision xel(2,nnodes_P2_0)
   complex(8) E_field_el(3,nnodes_P2_0)
   complex(8) H_field_el(3,nnodes_P2_0)
   complex(8) Ez_field_el_P3(nnodes_P3_0)
   double precision p2_p2(nnodes_P2_0,nnodes_P2_0)
   integer(8) i, j, j1
   integer(8) iel, ival
   integer(8) itrial, jtest, ui
   complex(8) vec_i(3), vec_j(3)
   complex(8) z_tmp1

   double precision mat_B(2,2), mat_T(2,2), det_b
!
!
!f2py intent(in) k_0, nval, nel, npt
!f2py intent(in) nnodes_P2, m_elnd_to_mshpt
!f2py intent(in) x, betas, soln_k1, type_el
!
!f2py depend(m_elnd_to_mshpt) nnodes_P2, nel
!f2py depend(x) npt
!f2py depend(betas) nval
!f2py depend(soln_k1) nnodes_P2, nval, nel
!f2py depend(type_el) nel
!
!f2py intent(out) overlap
!
!
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
!
   ui = stdout

!
   if ( nnodes_P2 .ne. 6 ) then
      write(ui,*) "EM_mode_en_int_v2: problem nnodes = ", nnodes_P2
      write(ui,*) "EM_mode_en_int_v2: nnodes should be equal to 6 !"
      write(ui,*) "EM_mode_en_int_v2: Aborting..."
      stop
   endif
!
   do ival=1,nval
      overlap1 = 0.0d0
      beta1 = betas(ival)
      do iel=1,nel
         typ_e = type_el(iel)
         if (typ_e .gt. 1) then
            do j=1,nnodes_P2
               j1 = m_elnd_to_mshpt(j,iel)
               nod_el_p(j) = j1
               xel(1,j) = x(1,j1)
               xel(2,j) = x(2,j1)
            enddo
!ccccccc
!       The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
!       maps the current triangle to the reference triangle.
            do i=1,2
               do j=1,2
                  mat_B(j,i) = xel(j,i+1) - xel(j,1)
               enddo
            enddo
            det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
            if (abs(det_b) .le. 1.0d-22) then
!c        if (abs(det_b) .le. 1.0d-8) then
               write(*,*) '?? EM_mode_energy_int_v2: Deter. = 0 :', det_b
               write(*,*) "xel = ", xel
               write(*,*) 'Aborting...'
               stop
            endif
!       We also need, is the matrix mat_T of the reverse transmation
!                (from reference to current triangle):
!       mat_T = inverse matrix of de mat_B
            mat_T(1,1) =  mat_B(2,2) / det_b
            mat_T(2,2) =  mat_B(1,1) / det_b
            mat_T(1,2) = -mat_B(1,2) / det_b
            mat_T(2,1) = -mat_B(2,1) / det_b
!       Note that if grad_i_0 is the gradient on the reference triangle,
!       then the gradient on the actual triangle is:
!       grad_i  = Transpose(mat_T)*grad_i0
!ccccccc
            do i=1,nnodes_P2
!         The components (E_x,E_y) of the mode ival
               do j=1,2
                  z_tmp1 = soln_k1(j,i,ival,iel)
                  E_field_el(j,i) = z_tmp1
               enddo
!         The component E_z of the mode ival. The FEM code uses the scaling:
!         E_z = C_IM_ONE* beta1 * \hat{E}_z
               j=3
               z_tmp1 = soln_k1(j,i,ival,iel)
               E_field_el(j,i) = z_tmp1 * C_IM_ONE* beta1
            enddo
!       E_z-field:
            do i=1,3
!         The longitudinal component at the vertices (P3 elements)
               j=3
               z_tmp1 = soln_k1(j,i,ival,iel)
               Ez_field_el_P3(i) = z_tmp1 * C_IM_ONE* beta1
            enddo
            do i=4,nnodes_P3_0
!         The longitudinal component at the edge nodes and interior node (P3 elements)
               j=3
               z_tmp1 = soln_k1(j,i+nnodes_P2-3,ival,iel)
               Ez_field_el_P3(i) = z_tmp1 * C_IM_ONE* beta1
            enddo
!
            call get_H_field_p3(k_0, beta1, mat_T,&
            &E_field_el, Ez_field_el_P3, H_field_el)
!
!       The matrix p2_p2 contains the overlap integrals between the P2-polynomial basis functions
            call find_overlaps_p2_p2(p2_p2, det_b)
            do itrial=1,nnodes_P2_0
               do i=1,3
                  z_tmp1 = E_field_el(i,itrial)
                  vec_i(i) = z_tmp1
               enddo
               do jtest=1,nnodes_P2_0
                  do i=1,3
                     z_tmp1 = H_field_el(i,jtest)
                     vec_j(i) = z_tmp1
                  enddo
!           Cross-product Z.(E^* X H) of E^*=vec_i and H=vec_j
                  z_tmp1 = vec_i(1) * vec_j(2) - vec_i(2) * vec_j(1)
                  overlap1 = overlap1 + z_tmp1 * p2_p2(itrial, jtest)
               enddo
            enddo
         endif
      enddo
      overlap(ival) = overlap1
   enddo
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
end subroutine EM_mode_energy_int_v2_wg
