
 !!!!!!!!!!!!!!!!


 !!!!!!!!!!!!!!!!

subroutine mat_el_v3 (xel, beta, c_tensor_el, rho_el,&
mat_K, mat_M)

   !!!!!!!!!!!!!!!!

   use numbatmod
   double precision xel(2,6)
   complex(8) beta
   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el

   !  Local variables

   double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
   complex(8) z_tmp1, z_tmp2
   complex(8) z_mat_xyz(6,6,3,3), z_beta, z_tensor
   integer(8) i, j, i_p, j_p, i_xyz,  j_xyz
   integer(8) i_u_xyz,  j_u_xyz, i_ind, j_ind
   integer(8) S_index(3,3)
   integer(8) debug

   !  Compute the Affine mappings from the current triangle to the
   !  reference unit triangle.
   !  Integration will be performed on the reference unit triangle


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c

   debug = 0



   !  S_index(i_xyz,i_u_xyz): index of S_{i_xyz,i_u_xyz} in the Voigt notation
   !  i_xyz represents the x,y, or z-derivative
   !  i_u_xyz represents the x,y, or z-field component

   !  S_xx => 1
   S_index(1,1) = 1
   !  S_yx => 6
   S_index(2,1) = 6
   !  S_zx => 5
   S_index(3,1) = 5

   !  S_xy => 6
   S_index(1,2) = 6
   !  S_yy => 2
   S_index(2,2) = 2
   !  S_zy => 4
   S_index(3,2) = 4

   !  S_xz => 5
   S_index(1,3) = 5
   !  S_yz => 4
   S_index(2,3) = 4
   !  S_zz => 3
   S_index(3,3) = 3

   do i=1,2
      do j=1,2
         mat_B(j,i) = xel(j,i+1) - xel(j,1)
      enddo
   enddo
   det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
   !  TEMPORARY CHANGE
   if (abs(det_b) .le. 1.0d-22) then
      !c      if (abs(det_b) .le. 1.0d-8) then
      write(*,*) '?? mat_el_v3: Determinant = 0 :', det_b
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

   do i=1,6
      do j=1,6
         do i_xyz=1,3
            do j_xyz=1,3
               z_mat_xyz(i,j,i_xyz,j_xyz) = 0
            enddo
         enddo
      enddo
   enddo

   !  z_mat_xyz: contains the overlap integrals of the x,y and z-derivatives
   do i=1,6
      do j=1,6
         do i_xyz=1,3
            do j_xyz=1,3
               z_mat_xyz(i,j,i_xyz,j_xyz) = 0
               if (i_xyz == 1 .and. j_xyz == 1) then
                  z_tmp1 = p2x_p2x(i,j)
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 1 .and. j_xyz == 2) then
                  z_tmp1 = p2x_p2y(i,j)
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 1 .and. j_xyz == 3) then
                  z_beta = C_IM_ONE* beta
                  z_tmp1 = p2_p2x(j,i) * z_beta
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 2 .and. j_xyz == 1) then
                  z_tmp1 = p2x_p2y(j,i)
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 2 .and. j_xyz == 2) then
                  z_tmp1 = p2y_p2y(i,j)
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 2 .and. j_xyz == 3) then
                  z_beta = C_IM_ONE* beta
                  z_tmp1 = p2_p2y(j,i) * z_beta
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 3 .and. j_xyz == 1) then
                  !  Conjugate
                  z_beta = - C_IM_ONE* beta
                  z_tmp1 = p2_p2x(i,j) * z_beta
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 3 .and. j_xyz == 2) then
                  !  Conjugate
                  z_beta = - C_IM_ONE* beta
                  z_tmp1 = p2_p2y(i,j) * z_beta
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               elseif (i_xyz == 3 .and. j_xyz == 3) then
                  z_beta = beta**2
                  z_tmp1 = p2_p2(i,j) * z_beta
                  z_mat_xyz(i,j,i_xyz,j_xyz) = z_tmp1
               endif
            enddo
         enddo
      enddo
   enddo

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
      !  Components of the displacement vector
      do i_u_xyz=1,3
         i_p = 3*(i-1) + i_u_xyz
         do j=1,6
            !  Components of the displacement vector
            do j_u_xyz=1,3
               j_p = 3*(j-1) + j_u_xyz
               !  Derivatives
               do i_xyz=1,3
                  i_ind = S_index(i_xyz,i_u_xyz)
                  !  Derivatives
                  do j_xyz=1,3
                     j_ind = S_index(j_xyz,j_u_xyz)
                     z_tensor = c_tensor_el(i_ind,j_ind)
                     z_tmp1 = z_mat_xyz(i,j,i_xyz,j_xyz)
                     if (i_u_xyz == 3) then
                        z_tmp1 = -C_IM_ONE* z_tmp1
                     endif
                     if (j_u_xyz == 3) then
                        z_tmp1 = C_IM_ONE* z_tmp1
                     endif
                     z_tmp2 = z_tmp1 * z_tensor
                     mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp2
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
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
   if (debug .eq. 1) then
      open(4,file="Output/xel.txt",status='unknown')
      do i=1,6
         write(4,"(I6,2(e20.10))") i, xel(1,i), xel(2,i)
      enddo
      close(4)
   endif

   return
end

