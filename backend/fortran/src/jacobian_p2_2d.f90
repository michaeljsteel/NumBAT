
!    Matrix B: Compute the Affine mappings from the current triangle to the
!     reference unit triangle N. Integration will be performed on that
!     triangle by Gaussian quadrature.

! This is for the isoparametric elements (never happens),
! where we use functions to construct transformation matrix
!
!    x_act = B X + (x_0, y_0)^t
!
! mat_B_0 is linear mapping from reference to actual
! mat_T is inverse mapping
!
subroutine jacobian_p2_2d (x_el, nnodes, p2_list,&
grad_p2_mat, x_act, det_jacobian, mat_B_0, mat_T, errco, emsg)
!
   use numbatmod
   integer(8) nnodes
   double precision x_el(2,nnodes)
   double precision p2_list(nnodes), grad_p2_mat(2,nnodes)
   double precision mat_B(2,2), mat_T(2,2), mat_B_0(2,2)
   double precision det_jacobian, x_act(2)
   integer(8) inode, i, j
   double precision phi, grad_phi
   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


!     32-but integers for BLAS and LAPACK
   integer(4) INFO_32, LDB_32, LDT_32
   integer(4) IPIV_32(2), NRHS_32, N_32

!
   ! do i=1,2
   !    x_act(i) = 0.0d0
   !    do j=1,2
   !       mat_B(i,j) = 0.0d0
   !    enddo
   ! enddo
   x_act = D_ZERO
   mat_B = D_ZERO

   do inode=1,nnodes
      do i=1,2
         phi = p2_list(inode)
         x_act(i) = x_act(i) + x_el(i,inode)*phi
         do j=1,2
            grad_phi = grad_p2_mat(j,inode)
            mat_B(i,j) = mat_B(i,j) + x_el(i,inode)*grad_phi
         enddo
      enddo
   enddo

   mat_B_0 = mat_B

   ! do i=1,2
   !    do j=1,2
   !       mat_B_0(j,i) = mat_B(j,i)
   !    enddo
   ! enddo


!       !  The order of the matrix mat_B
   N_32 = 2
!       !  The number of right hand sides
   NRHS_32 = N_32
!       !  The leading dimension of the array mat_B
   LDB_32 = N_32
!       !  The leading dimension of the array mat_T
   LDT_32 = N_32
!
! !     Initialisation for DGESV: mat_T = identity
!    do i=1,2
!       do j=1,2
!          mat_T(i,j) = 0.0d0
!       enddo
!       mat_T(i,i) = 1.0d0
!    enddo

   !  Initialisation for DGESV: mat_T = identity
   mat_T =  reshape( (/ 1.d0, 0.d0, &
                        0.d0, 1.d0  /), shape(mat_T), order=(/2,1/))


   call DGESV( N_32, NRHS_32, mat_B, LDB_32, IPIV_32,&
   &mat_T, LDT_32, INFO_32 )

     if(INFO_32 .ne. 0) then
      errco = NBERR_BAD_JACOBIAN
      write(emsg,*) 'Bad jacobian in jacobian_p1_2d'
      return
   endif

!     The value determinant can be obtained from the factorization P*L*U
   det_jacobian = 1
   do i=1,2
      if( (IPIV_32(i)-i) .eq. 0) then
         det_jacobian = det_jacobian*mat_B(i,i)
      else
         det_jacobian = -det_jacobian*mat_B(i,i)
      endif
   enddo

   if(abs(det_jacobian) .lt. 1.0d-20) then
      errco = NBERR_BAD_JACOBIAN
      write(emsg,*) 'Bad determinant in jacobian_p1_2d'
      do i=1,nnodes
         write(emsg,*) "i, x, y, x, = ", i, (x_el(j,i),j=1,2)
      enddo
      return

   endif

end
