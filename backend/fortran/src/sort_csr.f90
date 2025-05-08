! sort the elements within each row of a CSR rep
! (actually used to sort elements within each column of a csc rep)
subroutine sort_csr (n_dof, nonz, max_col_len, row_ind, col_ptr)

   use numbatmod
   use alloc
   type (NBError) nberr

   integer(8) n_dof, nonz, max_col_len
   integer(8) col_ptr(n_dof+1), row_ind(nonz)

   ! -----------------------------------------------

   integer(8) row_start, row_end, col_len
   integer(8) i, j, k

   integer(8), dimension(:), allocatable :: arr, istack, indx

   ! -----------------------------------------------

   call integer_alloc_1d(arr, max_col_len, 'arr', nberr)
   call integer_alloc_1d(istack, max_col_len, 'arr', nberr)
   call integer_alloc_1d(indx, max_col_len, 'arr', nberr)



   do i=1,n_dof
      row_start = col_ptr(i)
      row_end = col_ptr(i+1) - 1
      col_len = row_end - row_start + 1

      do j=row_start,row_end
         k = j - row_start + 1
         arr(k) = row_ind(j)
      enddo

      call sort_int (col_len, arr, indx, istack)

      do j=row_start,row_end
         k = j - row_start + 1
      enddo

      do j=row_start,row_end
         k = j - row_start + 1
         row_ind(j) = arr(indx(k))
      enddo

   enddo


   return
end

