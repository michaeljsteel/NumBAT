
!	Construction of the matrix of power flow
!	(integral of the z-component of the acoustic Poynting vector)
!

      subroutine mat_el_energy_rho (xel, rho_el, mat_P)


         implicit none
         double precision xel(2,6)
         complex(8) mat_P(6,6)
         complex(8) rho_el


         double precision p2_p2(6,6)
         double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
         double precision det_b
         complex(8)  z_tmp1
         integer(8) i, j
         integer(8) debug

         !  Compute the Affine mappings from the current triangle to the
         !   reference unit triangle.
         !  Integration will be performed on the reference unit triangle

         debug = 0


         do i=1,2
            do j=1,2
               mat_B(j,i) = xel(j,i+1) - xel(j,1)
            enddo
         enddo

         det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
         if (abs(det_b) .le. 1.0d-22) then
            !c      if (abs(det_b) .le. 1.0d-8) then
            write(*,*) '?? mat_el_energy_rho: Determinant = 0 :', det_b
            write(*,*) "xel = ", xel
            write(*,*) 'Aborting...'
            stop
         endif

         !   mat_T = Inverse of mat_B
         mat_T(1,1) = mat_B(2,2) / det_b
         mat_T(2,2) = mat_B(1,1) / det_b
         mat_T(1,2) = -mat_B(1,2) / det_b
         mat_T(2,1) = -mat_B(2,1) / det_b

         !	mat_T_tr = Tanspose(mat_T)

         mat_T_tr(1,1) = mat_T(1,1)
         mat_T_tr(1,2) = mat_T(2,1)
         mat_T_tr(2,1) = mat_T(1,2)
         mat_T_tr(2,2) = mat_T(2,2)

         call find_overlaps_p2_p2(p2_p2, det_b)

         do i=1,6
            do j=1,6
               mat_P(j,i) = 0
            enddo
         enddo

         !=================  Construction of the matrix of power flow   ==================
         !                 (integral of the z-component of the acoustic Poynting vector)
         do i=1,6
            do j=1,6
               z_tmp1 = p2_p2(i,j) * rho_el
               mat_P(i,j) = mat_P(i,j) + z_tmp1
            enddo
         enddo


         return
      end

