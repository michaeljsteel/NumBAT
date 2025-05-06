#include "numbat_decl.h"

! This one is written in CSC format

subroutine csc_make_col_ptr_loose (nel, n_entty, ety_tags, n_dof,  m_eqs, v_col_ptr, nonz_max)

   use numbatmod

   integer(8) nel, n_dof, n_entty, nonz_max
   integer(8) ety_tags (14,nel)
   integer(8) m_eqs(3,n_entty), v_col_ptr(n_dof+1)


   integer(8) k, k_el, active_dof, tag, i_nd, i_col, locdof
   integer(8) last_ct, this_ct


   v_col_ptr = 0


   !  Count references to each dof
   !  This is equivalent to counting how many elements a given entity falls on.

   !  Here v_col_ptr is just a convenient temporary memory holder.
   !  The contents is not related to its actual definition as the coloumn pointer

   do  k_el=1,nel           ! for every ety at every elt
      do i_nd=1,N_ENTITY_PER_EL
         tag = ety_tags(i_nd,k_el)      ! find its tag
         do locdof = 1,3                   ! for each of its 3 possible dof,
            active_dof = m_eqs(locdof, tag)    ! find the index of that dof, nonzero means active
            if (active_dof .ne. 0) v_col_ptr(active_dof) = v_col_ptr(active_dof)+1  ! count the number of times the dof is encountered
         enddo
      enddo
   enddo


   ! v_col_ptr now contains the number of mentions of each dof, ie the number of equations it is involved in
   ! This is more or less the number of triangles it lies on,
   ! The largest values occur at elt entities 5,6 or probably 7 which are vertices when many triangles meet
   ! This is usually of order 6-10


   !  Compressed Column Storage (CSC): determine the column pointer
   !   Each v_col_ptr(j)  = 1 + number of nonzeros left of column j
   !                      = v_col_ptr(j-1) + number of nonzeros column (j-1)
   !                                  where we set the virtual zero column as v_col_ptr(0)=1
   last_ct = v_col_ptr(1)
   v_col_ptr(1) = 1          !  The first column begins with element 1 (every dof has some self energy so this is always true )
   do i_col=2,n_dof+1
      this_ct = v_col_ptr(i_col)    ! # of neigbours of entity i_col (including self)

      ! Set v_col_ptr(i_col) by counting the possible elements in the previous column
      ! An upper bound for how many dof this one might interact with is:
      !     3 dof for each N_ENTITY_PER_EL entities in its own elt, + 3 dof for each of the (N_ENTITY_PER_EL-1) neighbour entities
      ! Some of these will be redundant by double counting, for example the edge and vertex nodes on two adjacent triangles
      !   which are themselves touching
      ! This is still _much_ smaller than interacting with every dof in the mesh
      v_col_ptr(i_col) = v_col_ptr(i_col-1) + 3*N_ENTITY_PER_EL + 3*(N_ENTITY_PER_EL-1)*(last_ct-1)
      last_ct = this_ct
   enddo

   nonz_max = v_col_ptr(n_dof+1) - 1

end subroutine


! Do this by pasing to a CSR routine with arrays in reverse order
! The matrix is symmetric so this is ok

! subroutine csc_make_col_ptr_tight (nel, n_entty, ety_tags, n_dof,  m_eqs, v_col_ptr, nonz_max)


! end subroutine



! This one is written in CSR format
! row/col names seem backward
! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine csr_length (n_msh_elts, n_entty, n_dof,  ety_tags, m_eqs, &
   col_ind, row_ptr, &  ! these names are swtiched from the call, but matched to the weird reverse naming in this file
   nonz_max, nonz, max_row_len, nberr)

   use numbatmod
   use alloc
   type(NBError) nberr

   integer(8) n_msh_elts, n_entty, n_dof
   integer(8) ety_tags(14,n_msh_elts)
   integer(8) m_eqs(3,n_entty)


   integer(8), dimension(:), allocatable, intent(inout) :: col_ind
   integer(8), dimension(:) :: row_ptr(n_dof+1)

   integer(8) nonz_max, nonz, max_row_len

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! --------------------------------------------

   integer(8), dimension(:), allocatable :: col_ind_0

   integer(8) i, j, j_nd, k, k1,  i_nd, i_tag, j_tag, i_locdof, j_locdof
   integer(8) k_el, i_dof, j_dof, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2, ui, stored
   integer(8) ct


   ui = stdout


   call integer_alloc_1d(col_ind_0, nonz_max, 'col_ind_0', nberr); RET_ON_NBERR(nberr)

   col_ind_0 = 0

   ! This code is expressed in CSR format
   ! The row shifts are known (which are actually those in v_col_ptr outside this function)
   ! Now we find the column index for every nonzero entry
   !  Determination of the column indices

   write(*,*) 'Total dof is ', n_entty, n_dof, nonz_max
   nonz = 0
   do k_el=1,n_msh_elts                    ! for each element

      do i_nd=1,N_ENTITY_PER_EL               !   and its 14 entities
         i_tag = ety_tags(i_nd, k_el)

         do i_locdof=1,3                     !   and their 3 potential dof
            i_dof = m_eqs(i_locdof, i_tag)   !   When nonzero, this is the row number for this dof

            if (i_dof .eq. 0) cycle     ! an inactive dof for this entity, go around again



            row_start = row_ptr(i_dof)          ! range of elt indices which are in this row
            row_end = row_ptr(i_dof+1) - 1

            !write(*,*) 'looking for partners of', k_el, i_nd, i_tag, i_locdof, &
             !  i_dof, 'in cols', row_start, row_end, row_end-row_start+1

            do j_nd=1,N_ENTITY_PER_EL
               j_tag = ety_tags(j_nd,k_el)

               do j_locdof=1,3
                  j_dof = m_eqs(j_locdof, j_tag)

                  if (j_dof .eq. 0) cycle

                  !  Store the entry (i_dof,j_dof) if it's not already found
                  stored = 0

                  do k=row_start,row_end
                     if (col_ind_0(k) .eq. 0) then  !goto 20 ! if we find a zero, we've seen all of the nonzeros, and none were a match
                        !  No entry exists for (i_dof,j_dof); create new one
                        nonz = nonz + 1

                        col_ind_0(k) = j_dof
                        stored = 1
                        exit
                     endif

                     if (col_ind_0(k) .eq. j_dof) then !goto 30  ! already stored, bail out
                        stored=1
                        exit
                     endif
                  enddo

                  if (stored .eq. 0) then ! shouldn't have got here
                     emsg = "csr_length: There is a problem with row/col indexing!"
                     call nberr%set(NBERROR_188, emsg)
                     return
                  endif

               enddo
            enddo
         enddo
      enddo
   enddo

   write(*,*) 'csrlen 2'


   ! The full coupling matrix is dimension n_dof x n_dof
   ! At most nonz_max of these can be nonzero based on the pairs lying on adjacent elements
   ! We now have a valid CSR indexing of this subset of nonz_max pairs

   ! But more of these can be eliminated. Not quite sure why.


   !  squeeze away the zero entries
   !  added so as to handle more type of domains/meshes

!    do i=1,n_dof-1
!       row_start = row_ptr(i)
!       row_end = row_ptr(i+1) - 1

!       do j=row_start,row_end
!          if(col_ind_0(j) .eq. 0) then   ! if one zero, then all the rest of this row will be zeros, so remove them all
!             row_start2 = row_ptr(i) + j - row_start
!             row_ptr(i+1) = row_start2         ! bring the start of the next row forward
!             row_end2 = row_ptr(i+2) - 1       ! find the end of that next row
!             do k=row_end+1,row_end2
!                k1 = row_start2 + k - (row_end+1)  ! shuufle the columns in that row forward into the empty space
!                col_ind_0(k1) = col_ind_0(k)
!                col_ind_0(k) = 0
!             enddo
!             goto 40
!          endif
!       enddo

! 40    continue
!    enddo


   write(*,*) 'csrlen 3'
   do i=1,n_dof-1
      row_start = row_ptr(i)
      row_end = row_ptr(i+1) - 1

      do j=row_start,row_end
         if(col_ind_0(j) .ne. 0) cycle   ! this elt is busy, check the next

         ! we've found a zero element in this row, all the rest of this row will be zeros, so remove them all

         row_start2 = row_ptr(i) + j - row_start
         row_ptr(i+1) = row_start2         ! bring the start of the next row forward
         row_end2 = row_ptr(i+2) - 1       ! find the end of that next row
         do k=row_end+1,row_end2
            k1 = row_start2 + k - (row_end+1)  ! shuufle the columns in that row forward into the empty space
            col_ind_0(k1) = col_ind_0(k)
            col_ind_0(k) = 0
         enddo

         exit  ! this row is done, go back to the next

      enddo


   enddo



   ! squeeze the last row
   i = n_dof
   row_start = row_ptr(i)
   row_end = row_ptr(i+1) - 1
   do j=row_start,row_end
      if(col_ind_0(j) .eq. 0) then
         row_start2 = row_ptr(i) + j - row_start
         row_ptr(i+1) = row_start2
         exit
      endif
   enddo


   ! Find the longest row
   max_row_len = 0
   do i=1,n_dof
      row_start = row_ptr(i)
      row_end = row_ptr(i+1) - 1
      row_len = row_end - row_start + 1
      if (row_len .gt. max_row_len) max_row_len = row_len
   enddo

   write(*,*) 'csrlen 5', nonz


   ! Now we know nonz
   call integer_alloc_1d(col_ind, nonz, 'col_ind', nberr); RET_ON_NBERR(nberr)

   write(*,*) 'csrlen 6'

   col_ind(1:nonz) = col_ind_0(1:nonz)

   write(*,*) 'csrlen 7'

!   deallocate(col_ind_0)

end
