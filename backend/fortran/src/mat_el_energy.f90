

 !	Construction of the matrix of power flow
 !	(integral of the z-component of the acoustic Poynting vector)
 !

subroutine mat_el_energy (xel, beta, c_tensor_el, mat_P)

   use numbatmod
   double precision xel(2,6)
   complex(8) beta
   complex(8) mat_P(18,18)
   complex(8) c_tensor_el(6,6)

   !     Local variables

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
   complex(8) z_tmp1
   integer(8) i, j, i_p, j_p, i_xyz,  j_xyz
   integer(8) debug

   !    Compute the Affine mappings from the current triangle to the
   !     reference unit triangle.
   !    Integration will be performed on the reference unit triangle



   debug = 0


   do i=1,2
      do j=1,2
         mat_B(j,i) = xel(j,i+1) - xel(j,1)
      enddo
   enddo
   det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
   !       !  TEMPORARY CHANGE
   if (abs(det_b) .le. 1.0d-22) then
      !c      if (abs(det_b) .le. 1.0d-8) then
      write(*,*) '?? mat_el_energy: Determinant = 0 :', det_b
      write(*,*) "xel = ", xel
      write(*,*) 'Aborting...'
      stop
   endif
   !     mat_T = Inverse of mat_B
   mat_T(1,1) = mat_B(2,2) / det_b
   mat_T(2,2) = mat_B(1,1) / det_b
   mat_T(1,2) = -mat_B(1,2) / det_b
   mat_T(2,1) = -mat_B(2,1) / det_b
   !
   !
   !	mat_T_tr = Tanspose(mat_T)
   !
   mat_T_tr(1,1) = mat_T(1,1)
   mat_T_tr(1,2) = mat_T(2,1)
   mat_T_tr(2,1) = mat_T(1,2)
   mat_T_tr(2,2) = mat_T(2,2)

   call mat_p2_p2(p2_p2, det_b)
   call mat_p2_p2x (p2_p2x, mat_T_tr, det_b)
   call mat_p2_p2y (p2_p2y, mat_T_tr, det_b)
   !
   do i=1,18
      do j=1,18
         mat_P(j,i) = 0
      enddo
   enddo

   !=================  Construction of the matrix of power flow   ==================
   !                   (integral of the z-component of the acoustic Poynting vector)

   do i=1,6
      do i_xyz=1,3
         i_p = 3*(i-1) + i_xyz
         do j=1,6
            do j_xyz=1,3
               j_p = 3*(j-1) + j_xyz
               if (i_xyz == 1 .and. j_xyz == 1) then
                  z_tmp1 = p2_p2(i,j) * beta * C_IM_ONE
                  ! C(5,5) * u_x * S_zx
                  z_tmp1 = z_tmp1 * c_tensor_el(5,5)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 1 .and. j_xyz == 3) then
                  z_tmp1 = p2_p2x(i,j)
                  !   C(5,5) * u_x * S_xz
                  z_tmp1 = z_tmp1 * c_tensor_el(5,5)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 2 .and. j_xyz == 2) then
                  z_tmp1 = p2_p2(i,j) * beta * C_IM_ONE
                  !   C(4,4) * u_x * S_zy
                  z_tmp1 = z_tmp1 * c_tensor_el(4,4)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 2 .and. j_xyz == 3) then
                  z_tmp1 = p2_p2y(i,j)
                  !  C(4,4) * u_x * S_yz
                  z_tmp1 = z_tmp1 * c_tensor_el(4,4)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 3 .and. j_xyz == 1) then
                  z_tmp1 = p2_p2x(i,j)
                  !  C(3,1) * u_x * S_xx
                  z_tmp1 = z_tmp1 * c_tensor_el(3,1)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 3 .and. j_xyz == 2) then
                  z_tmp1 = p2_p2y(i,j)
                  !  C(3,2) * u_x * S_yy
                  z_tmp1 = z_tmp1 * c_tensor_el(3,2)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1

               elseif (i_xyz == 3 .and. j_xyz == 3) then
                  z_tmp1 = p2_p2(i,j) * beta * C_IM_ONE
                  !  C(3,3) * u_x * S_zz
                  z_tmp1 = z_tmp1 * c_tensor_el(3,3)
                  mat_P(i_p,j_p) = mat_P(i_p,j_p) - z_tmp1
               endif
            enddo
         enddo
      enddo
   enddo



   if (debug .eq. 1) then
      open(4,file="Output/mat_P.txt",status='unknown')
      do i=1,18
         do j=1,18
            write(4,"(2(I6),6(e20.10))") i,j, mat_P(j,i), mat_P(i,j), mat_P(i,j) - conjg(mat_P(j,i))
         enddo
      enddo
      close(4)
   endif

   if (debug .eq. 1) then
      open(4,file="Output/c_tensor_el.txt",status="unknown")
      do j=1,6
         do i=1,6
            write(4,"(2(I6),2(e20.10))") i,j, c_tensor_el(i,j)
         enddo
      enddo
      close(4)
   endif


   if (debug .eq. 1) then
      open(4,file="Output/xel.txt",status='unknown')
      do i=1,6
         write(4,"(I6,2(e20.10))") i, xel(1,i), xel(2,i)
      enddo
      close(4)
   endif


   return
end

