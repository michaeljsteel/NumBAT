

subroutine mat_el_v2 (xel, beta, c_tensor_el, rho_el,&
mat_K, mat_M)


   implicit none
   double precision xel(2,6)
   complex(8) beta
   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el

   !  Local variables

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
   complex(8) z_tmp1
   integer(8) i, j, i_p, j_p, i_xyz,  j_xyz
   integer(8) debug

   !  Compute the Affine mappings from the current triangle to the
   !  reference unit triangle.
   !  Integration will be performed on the reference unit triangle


   !!!!!!!!!!!c

   debug = 0



   do i=1,2
      do j=1,2
         mat_B(j,i) = xel(j,i+1) - xel(j,1)
      enddo
   enddo

   det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
   !  TEMPORARY CHANGE
   if (abs(det_b) .le. 1.0d-22) then
      !c      if (abs(det_b) .le. 1.0d-8) then
      write(*,*) '?? mat_el_v2: Determinant = 0 :', det_b
      write(*,*) "xel = ", xel
      write(*,*) 'Aborting...'
      stop
   endif

   !  mat_T = Inverse of mat_B
   mat_T(1,1) = mat_B(2,2) / det_b
   mat_T(2,2) = mat_B(1,1) / det_b
   mat_T(1,2) = -mat_B(1,2) / det_b
   mat_T(2,1) = -mat_B(2,1) / det_b


   !  mat_T_tr = Tanspose(mat_T)

   mat_T_tr(1,1) = mat_T(1,1)
   mat_T_tr(1,2) = mat_T(2,1)
   mat_T_tr(2,1) = mat_T(1,2)
   mat_T_tr(2,2) = mat_T(2,2)

   call find_overlaps_p2_p2(p2_p2, det_b)
   call find_overlaps_p2_p2x (p2_p2x, mat_T_tr, det_b)
   call find_overlaps_p2_p2y (p2_p2y, mat_T_tr, det_b)
   call find_overlaps_p2x_p2x (p2x_p2x, mat_T_tr, det_b)
   call find_overlaps_p2x_p2y (p2x_p2y, mat_T_tr, det_b)
   call find_overlaps_p2y_p2y (p2y_p2y, mat_T_tr, det_b)

   do i=1,18
      do j=1,18
         mat_K(j,i) = 0
         mat_M(j,i) = 0
      enddo
   enddo

   !=================  Construction of the matrix mat_M =================
   !  Integral [rho * P(i) * P(i)]
   do i=1,6
      !  The components x, y and z
      do i_xyz=1,3
         i_p = 3*(i-1) + i_xyz
         do j=1,6
            j_xyz = i_xyz
            j_p = 3*(j-1) + j_xyz
            z_tmp1 = p2_p2(i,j) * rho_el
            mat_M(i_p,j_p) = mat_M(i_p,j_p) + z_tmp1
         enddo
      enddo
   enddo

   !=================  Construction of the matrix mat_K =================
   !  Integral [K_{ij} = Gradient_s(conjg(P_k(i))) x c_tensor x Gradient_s(P_k(j))], where k=x,y,z
   !  Reference: see Eqs. (7) and (8) in:
   !  A.-C. Hladky-Hennion
   !  "Finite element analysis of the propagation of acoustic waves in waveguides,"
   !  Journal of Sound and Vibration, vol. 194, no. 2, pp. 119-136, 1996.
   do i=1,6
      do i_xyz=1,3
         i_p = 3*(i-1) + i_xyz
         do j=1,6
            j_xyz = i_xyz
            j_p = 3*(j-1) + j_xyz
            if (i_xyz == 1) then
               !  Overlap: row 1 of [C]*[B] ###########
               z_tmp1 = p2x_p2x(i,j)
               !  C(1,1) * S_xx * conjg(S_xx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 5 of [C]*[B]  ###########
               z_tmp1 = p2_p2(i,j) * beta**2
               !  C(5,5) * S_zx * conjg(S_zx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 6 of [C]*[B]  ###########
               z_tmp1 = p2y_p2y(i,j)
               !  C(6,6) * S_yx * conjg(S_yx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,6)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
               !  Overlap: row 2 of [C]*[B]  ###########
               z_tmp1 = p2y_p2y(i,j)
               !  C(2,2) * S_yy * conjg(S_yy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 4 of [C]*[B]  ###########
               z_tmp1 = p2_p2(i,j) * beta**2
               !  C(4,4) * S_zy * conjg(S_zy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 6 of [C]*[B]  ###########
               z_tmp1 = p2x_p2x(i,j)
               !  C(6,6) * S_xy * conjg(S_xy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
               !  Overlap: column 3 of [C]*[B]  ###########
               z_tmp1 = p2_p2(i,j) * beta**2
               !  C(3,3) * S_zz * conjg(S_zz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 4 of [C]*[B]  ###########
               z_tmp1 = p2y_p2y(i,j)
               !  C(4,4) * S_yz * conjg(S_yz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,4)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 5 of [C]*[B]  ###########
               z_tmp1 = p2x_p2x(i,j)
               !  C(5,5) * S_xz * conjg(S_xz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
         enddo
      enddo
   enddo

   do i=1,6
      do i_xyz=1,3
         i_p = 3*(i-1) + i_xyz
         do j=1,6
            if (i_xyz == 1) then
               !  Overlap: row 1 of [C]*[B]  ###########
               j_xyz = 2
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2x_p2y(i,j)
               !  C(1,2) * S_xx * conjg(S_yy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 6 of [C]*[B]  ###########
               j_xyz = 2
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2x_p2y(j,i)
               !  C(6,6) * S_yx * conjg(S_xy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,6)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
               !  Overlap: row 2 of [C]*[B]  ###########
               j_xyz = 3
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2y(j,i) * (-beta)
               !  C(2,3) * S_yy * conjg(S_zz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 4 of [C]*[B]  ###########
               j_xyz = 3
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2y(i,j) * beta
               !  C(4,4) * S_zy * conjg(S_yz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
               !  Overlap: row 3 of [C]*[B]  ###########
               j_xyz = 1
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2x(i,j) * (-beta)
               !  C(3,1) * S_zz * conjg(S_xx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 5 of [C]*[B]  ###########
               j_xyz = 1
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2x(j,i) * beta
               !  C(5,5) * S_xz * conjg(S_zx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
         enddo
      enddo
   enddo

   do i=1,6
      do i_xyz=1,3
         i_p = 3*(i-1) + i_xyz
         do j=1,6
            if (i_xyz == 1) then
               !  Overlap: row 1 of [C]*[B]  ###########
               j_xyz = 3
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2x(j,i) * (-beta)
               !  C(1,3) * S_xx * conjg(S_zz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 5 of [C]*[B]  ###########
               j_xyz = 3
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2x(i,j) * beta
               !  C(5,5) * S_zx * conjg(S_xz)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
               !  Overlap: row 2 of [C]*[B]  ###########
               j_xyz = 1
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2x_p2y(j,i)
               !  C(2,1) * S_yy * conjg(S_xx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 6 of [C]*[B]  ###########
               j_xyz = 1
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2x_p2y(i,j)
               !  C(6,6) * S_xy * conjg(S_yx)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
               !  Overlap: row 3 of [C]*[B]  ###########
               j_xyz = 2
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2y(i,j) * (-beta)
               !  C(3,2) * S_zz * conjg(S_yy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
               !  Overlap: row 4 of [C]*[B]  ###########
               j_xyz = 2
               j_p = 3*(j-1) + j_xyz
               z_tmp1 = p2_p2y(j,i) * beta
               !  C(4,4) * S_yz * conjg(S_zy)
               z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,4)
               mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
         enddo
      enddo
   enddo
   !c
   !c
   if (debug .eq. 1) then
      open(4,file="Output/mat_K.txt",status='unknown')
      do i=1,18
         do j=1,18
            write(4,"(2(I6),6(e20.10))") i,j, mat_K(j,i),&
            &mat_K(i,j), mat_K(i,j) - conjg(mat_K(j,i))
         enddo
      enddo
      close(4)
   endif
   !c
   if (debug .eq. 1) then
      open(4,file="Output/mat_M.txt",status='unknown')
      do i=1,18
         do j=1,18
            write(4,"(2(I6),6(e20.10))") i,j, mat_M(j,i),&
            &mat_M(i,j), mat_M(i,j) - conjg(mat_M(j,i))
         enddo
      enddo
      close(4)
   endif
   !c
   if (debug .eq. 1) then
      open(4,file="Output/rho_el.txt",status="unknown")
      write(4,"(2(e20.10))") rho_el
      close(4)
      open(4,file="Output/c_tensor_el.txt",status="unknown")
      do j=1,6
         do i=1,6
            write(4,"(2(I6),2(e20.10))") i,j, c_tensor_el(i,j)
         enddo
      enddo
      close(4)
   endif
   !c
   !c
   if (debug .eq. 1) then
      open(4,file="Output/xel.txt",status='unknown')
      do i=1,6
         write(4,"(I6,2(e20.10))") i, xel(1,i), xel(2,i)
      enddo
      close(4)
   endif

   !c      write(*,*)
   !c      write(*,*) "mat_el_v2: Aborting..."
   !c      stop


   return
end

