! TODO: ex_arr, ex_istack are unused
!       remove from elastic call and then drop them


subroutine sort_csr (n_dof, nonz, max_row_len, col_ind, row_ptr, indx)

   use numbatmod
   use alloc
   type (NBError) nberr

   integer(8) n_dof, nonz, max_row_len
   integer(8) row_ptr(n_dof+1), col_ind(nonz)
   integer(8) indx(max_row_len)

   ! -----------------------------------------------

   integer(8) row_start, row_end, row_len
   integer(8) i, j, k

   integer(8), dimension(:), allocatable :: arr, istack

   ! -----------------------------------------------

   call integer_alloc_1d(arr, max_row_len, 'arr', nberr)
   call integer_alloc_1d(istack, max_row_len, 'arr', nberr)


   do i=1,n_dof
      row_start = row_ptr(i)
      row_end = row_ptr(i+1) - 1
      row_len = row_end - row_start + 1

      do j=row_start,row_end
         k = j - row_start + 1
         arr(k) = col_ind(j)
      enddo

      call sort_int (row_len, arr, indx, istack)

      do j=row_start,row_end
         k = j - row_start + 1
      enddo

      do j=row_start,row_end
         k = j - row_start + 1
         col_ind(j) = arr(indx(k))
      enddo

   enddo


   return
end

