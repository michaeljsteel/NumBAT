
! Find affine transformation from reference triangle to actual triangle
! mat_B_0 is linear mapping from reference to actual
! mat_T is inverse mapping

subroutine jacobian_p1_2d (x_ref, x_el, nnodes, x_act, det_jacobian, &
   mat_B_0, mat_T, errco, emsg)

   use numbatmod

   integer(8) nnodes
   double precision x_ref(2)                   ! coordinate in the reference triangle
   double precision x_act(2)                   ! coordinate in the actual triangle
   double precision x_el(2,nnodes)             ! coordinates of the element P2 nodes

   double precision mat_B(2,2), mat_T(2,2), mat_B_0(2,2)
   double precision det_jacobian

   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


   integer(8) i, j

   !  32-bit integers for BLAS and LAPACK
   integer(4) INFO_32, LDB_32, LDT_32
   integer(4) IPIV_32(2), NRHS_32, N_32

   ! Compute the Affine mappings from the current Tetrahedron to the
   !  reference unit Tetrahedron N. Integration will be performed on that
   !  Tetrahedron by Gaussian quadrature.

   ! x_act = B X + (x_0, y_0, z_0)^t

   ! B=[a2-a1, a3-a1], columns of the vectors between the P1 nodes
   do i=1,2
      !do j=1,2
      !mat_B(j,i) = x_el(j,i+1) - x_el(j,1)
      mat_B(:,i) = x_el(:,i+1) - x_el(:,1)
      !enddo
   enddo

   mat_B_0 = mat_B
   ! do i=1,2
   ! do j=1,2
   !    mat_B_0(j,i) = mat_B(j,i)
   ! enddo
   ! enddo


   N_32 = 2
   NRHS_32 = N_32
   LDB_32 = N_32
   LDT_32 = N_32

   ! Warp by B and translate reference origin to node 1 of actual triangle
   !   x_act = B X + (x_0, y_0, z_0)^t

   !  Initialisation for DGEMV x_act = x_el(1:2,1)
   !do i=1,2
   !   x_act(i) = x_el(i,1)
   !enddo
   x_act = x_el(:,1)          ! (x_0, y_0, z_0)^t

   ! x_act = B x + x_act
   call DGEMV('No transpose', N_32, N_32, D_ONE, mat_B, LDB_32, x_ref, 1, D_ONE, x_act, 1)


   ! do i=1,2
   ! do j=1,2
   !    mat_T(i,j) = 0.0d0
   ! enddo
   ! mat_T(i,i) = 1.0d0
   ! enddo

!  Initialisation for DGESV: mat_T = identity
   mat_T =  reshape( (/ 1.d0, 0.d0, &
                        0.d0, 1.d0  /), shape(mat_T), order=(/2,1/))



   ! mat_T = inv(mat_B)
   call DGESV( N_32, NRHS_32, mat_B, LDB_32, IPIV_32, mat_T, LDT_32, INFO_32 )
!
   if(INFO_32 .ne. 0) then
      errco = NBERR_BAD_JACOBIAN
      write(emsg,*) 'Bad jacobian in jacobian_p1_2d'
      return
   endif

!  The value determinant can be obtained from the factorization P*L*U
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
