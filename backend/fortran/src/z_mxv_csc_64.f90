
! Compute the product of a matrix in CSC format by a vector
! vect2 = mat * vect1

!TODO: Is there not a standard library function in SparseSuite for this?

subroutine z_mxv_csc (n_dof, vect1, vect2, nonz, row_ind, col_ptr, mat)

   use numbatmod
   integer(8) n_dof, nonz
   integer(8) row_ind(nonz), col_ptr(n_dof+1)
   complex(8) mat(nonz)
   complex(8) vect1(n_dof), vect2(n_dof)

   integer(8) i, j, k, col_start, col_end, i_base

   vect2 = C_ZERO

!  valpr.f has changed the CSC indexing to 0-based indexing
!  so we must add 1 to the CSD row_pointer row_ind
   i_base = 1


   do i=1,n_dof  !  Column index
      col_start = col_ptr(i) + i_base
      col_end = col_ptr(i+1) - 1 + i_base
      do j=col_start,col_end
         k = row_ind(j) + i_base  !  Row number
         !  mat(j) = value of the matrix entry (k,i)
         vect2(k) = vect2(k) + mat(j)*vect1(i)
      enddo
   enddo

end
