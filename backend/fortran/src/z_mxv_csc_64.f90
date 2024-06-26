
!***********************************************************************
!
!    Compute the product of a matrix in CSC format by a vector
!
!***********************************************************************

!TODO: Is there not a stnadard library function in SparseSuite for this?

subroutine z_mxv_csc (neq, vect1, vect2, nonz, row_ind, col_ptr, mat)

   use numbatmod
   integer(8) neq, nonz
   integer(8) row_ind(nonz), col_ptr(neq+1)
   complex(8) mat(nonz)
   complex(8) vect1(neq), vect2(neq)

   integer(8) i, j, k, col_start, col_end, i_base

   ! do i=1,neq
   !    vect2(i) = 0.d0
   ! enddo
   vect2 = C_ZERO

!     valpr.f has changed the CSC indexing to 0-based indexing
!     so we must add 1 to the CSD row_pointer row_ind
   i_base = 1


   do i=1,neq  ! Column index
      col_start = col_ptr(i) + i_base
      col_end = col_ptr(i+1) - 1 + i_base
      do j=col_start,col_end
         k = row_ind(j) + i_base  ! Row number
         !  mat(j) = value of the matrix entry (k,i)
         vect2(k) = vect2(k) + mat(j)*vect1(i)
      enddo
   enddo

   return
end
