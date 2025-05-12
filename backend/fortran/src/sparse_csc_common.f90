
subroutine store_dof_in_csc_row_index(row_lo, row_hi, j_dof, row_ind_tmp, n_nonz, n_nonz_max, nberr)

   use numbatmod
   integer(8) row_lo, row_hi, j_dof, n_nonz, n_nonz_max
   type(NBError) nberr

   integer(8) :: row_ind_tmp(n_nonz_max)
   logical found
   character(len=EMSG_LENGTH) :: emsg
   integer(8) k

   found = .false.
   do k=row_lo,row_hi

      ! We've run out of nonzero entries in the row indices, and haven't found this element yet
      ! So create one
      if(row_ind_tmp(k) .eq. 0) then
         found = .true.
         n_nonz =n_nonz + 1

         if (n_nonz .le. n_nonz_max) then ! We have an empty slot. Go!
            row_ind_tmp(k) = j_dof
            exit
         else  ! should never happen if make_provisional is correct
            write(emsg, *) "csr_length_AC:n_nonz > n_nonz_max: ",&
               n_nonz .gt. n_nonz_max
            call nberr%set(NBERR_BAD_SPARSE_FINAL, emsg)
            return

         endif


      endif

      ! Entry already exists, nothing to do
      if(row_ind_tmp(k) .eq. j_dof) then
         found = .true.
         exit
      endif

   enddo

   if (found .neqv. .true.) then
      write(emsg, *) "csr length looking for missing spot"
      call nberr%set(NBERR_BAD_SPARSE_FINAL, emsg)
      return
   endif
end subroutine