! TODO: this seems to normalise by the wrong factor: 1/|z|^2 instead of 1/|z|
! Is this ever actually used?

subroutine normalise_fields (n_modes, nel, nnodes, soln_k1, soln_k2, mat_overlap)
   implicit none
   integer(8) n_modes, nel, nnodes
   complex(8) soln_k1(3,nnodes+7,n_modes,nel)
   complex(8) soln_k2(3,nnodes+7,n_modes,nel)
   complex(8) mat_overlap(n_modes,n_modes)

   !integer(8) i,j
   integer(8) ival
   complex(8) z_tmp1, z_tmp2

   double precision, parameter :: min_abs = 1.0d-8

   do ival=1,n_modes
      !do iel=1,nel
         z_tmp1 = sqrt(mat_overlap(ival,ival))

         ! Don't bother rescaling fields which are super tiny and therefore probably broken anyway
         if (abs(z_tmp1) .gt. min_abs) then
!              z_tmp2 =  1.0d0/z_tmp1
            z_tmp2 =  1.0d0/z_tmp1**2
      !      do i=1,nnodes+7
      !         do j=1,3
      !            soln_k1(j,i,ival,iel) = soln_k1(j,i,ival,iel)
      !            soln_k2(j,i,ival,iel) = soln_k2(j,i,ival,iel) * z_tmp2
       !        enddo
        !    enddo
         !   soln_k1(j,i,ival,iel) = soln_k1(j,i,ival,iel)
            soln_k2(:,:,ival,:) = soln_k2(:,:,ival,:) * z_tmp2


         endif

      !enddo
   enddo

   return
end subroutine

