#include "numbat_decl.h"

! This one is written in CSC format

subroutine csc_make_col_ptr_loose (nel, n_ddl, ety_tags, neq,  m_global_dofs, v_col_ptr, nonz_max)

   use numbatmod

   integer(8) nel, neq, n_ddl, nonz_max
   integer(8) ety_tags (14,nel)
   integer(8) m_global_dofs(3,n_ddl), v_col_ptr(neq+1)


   integer(8) k, iel, active_dof, tag, i_nd, i_eq, dof
   integer(8) last_ct, this_ct, maxmeq


   v_col_ptr = 0


   !  Count references to each dof
   !  Here v_col_ptr is just a convenient temporary memory holder.
   !  The contents is not related to its actual definition.

   do  iel=1,nel           ! for every ety at every elt
      do i_nd=1,NDDL_0_EM
         tag = ety_tags(i_nd,iel)      ! find its tag
         do dof = 1,3                   ! for each of its 3 possible dof,
            active_dof = m_global_dofs(dof, tag)    ! find the index of that dof, nonzero means active
            if (active_dof .ne. 0) v_col_ptr(active_dof) = v_col_ptr(active_dof)+1  ! count the number of times the dof is encountered
            if (v_col_ptr(active_dof) .eq. 8) write(*,*) 'An 8 eqn dof', iel, i_nd, tag, dof, active_dof
         enddo
      enddo
   enddo


   ! v_col_ptr now contains the number of mentions of each dof, ie the number of equations it is involved in
   ! This is more or less the number of triangles it lies on,
   ! The largest values occur at elt entities 5,6 or probably 7 which are vertices when many triangles meet
   ! This is usually of order 6-10

   write(*,*) 'vptrs a', (v_col_ptr(k), k=1,20)

   !  NO: Compressed Row Storage (CRS): determine the row pointer
   !  Compressed Column Storage (CSC): determine the column pointer
   !   Each v_col_ptr(j)  = v_col_ptr(j-1) + number of nonzeros in the new column
   last_ct = v_col_ptr(1)
   v_col_ptr(1) = 1          !  The first column begins with element 1 (every dof has some self energy so this is always true )
   do i_eq=2,neq+1
      this_ct = v_col_ptr(i_eq)

      !An upper bound for how many dof this one might interact with is:
      !     3 dof for each NDDL_0_EM entities in its own elt, + 3 dof for each of the (NDDL_0_EM-1) entities
      ! This is still _much_ smaller than interacting with every dof in the mesh
      v_col_ptr(i_eq) = v_col_ptr(i_eq-1) + 3*NDDL_0_EM + 3*(NDDL_0_EM-1)*(last_ct-1)
      last_ct = this_ct
   enddo

   nonz_max = v_col_ptr(neq+1) - 1

   write(*,*) 'vptrs b', (v_col_ptr(k), k=1,10)
end


! This one is written in CSR format
! row/col names seem backward
! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine csr_length (n_msh_elts, n_ddl, neq,  ety_tags, m_global_dofs, &
   col_ind, row_ptr, &  ! these names are swtiched from the call, but matched to the weird reverse naming in this file
   nonz_max, nonz, max_row_len, debug, errco, emsg)

   use numbatmod
   use alloc

   integer(8) n_msh_elts, n_ddl, neq
   integer(8) ety_tags(14,n_msh_elts)
   integer(8) m_global_dofs(3,n_ddl)


   integer(8), dimension(:), allocatable, intent(inout) :: col_ind
   integer(8), dimension(:) :: row_ptr(neq+1)

   integer(8) nonz_max, nonz, max_row_len

   integer(8) debug

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! --------------------------------------------

   integer(8), dimension(:), allocatable :: col_ind_0

   integer(8) i, j, j_nd, k, k1, dof, j_ddl, i_nd, i_tag, j_tag, t_dof, i_locdof, j_locdof
   integer(8) iel, i_dof, ip, j_dof, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2, ui


   ui = stdout


   call integer_alloc_1d(col_ind_0, nonz_max, 'col_ind_0', errco, emsg); RETONERROR(errco)

   col_ind_0 = 0

   ! This code is expressed in CSR format
   ! The row shifts are known (which are actually those in v_col_ptr outside this function)
   ! Now we find the column index for every nonzero entry
   !  Determination of the column indices

   nonz = 0
   do iel=1,n_msh_elts                    ! for each element

      do i_nd=1,NDDL_0_EM               !   and its 14 entities
         i_tag = ety_tags(i_nd, iel)

         do i_locdof=1,3                     !   and their 3 potential dof
            i_dof = m_global_dofs(i_locdof, i_tag)   !   When nonzero, this is the row number for this dof

            if (i_dof .eq. 0) cycle     ! an inactive dof for this entity, go around again

            row_start = row_ptr(i_dof)          ! range of elt indices which are in this row
            row_end = row_ptr(i_dof+1) - 1

            write(*,*) 'rows', iel, i_nd, i_tag, i_locdof, i_dof, row_start, row_end, row_end-row_start+1

            do j_nd=1,NDDL_0_EM
               j_tag = ety_tags(j_nd,iel)

               do j_locdof=1,3
                  j_dof = m_global_dofs(j_locdof, j_tag)

                  if (j_dof .eq. 0) cycle

                     !  Search if the entry (i_dof,j_dof) is already stored
                     do k=row_start,row_end
                        if (col_ind_0(k) .eq. 0) goto 20 ! if we find a zero, we've seen all of the nonzeros, and none were a match
                        if (col_ind_0(k) .eq. j_dof) goto 30
                     enddo

                     ! bail out
                     emsg = "csr_length: There is a problem with row/col indexing!"
                     errco = NBERROR_118
                     return

20                   continue
                     !  No entry exists for (i_dof,j_dof); create new one
                     nonz = nonz + 1

                     col_ind_0(k) = j_dof

30                   continue
                  !endif

               enddo
            enddo
            !endif
         enddo
      enddo
   enddo


!  squeeze away the zero entries
!  added so as to handle more type of domains/meshes

   do i=1,neq-1
      row_start = row_ptr(i)
      row_end = row_ptr(i+1) - 1
      do j=row_start,row_end
         if(col_ind_0(j) .eq. 0) then
            row_start2 = row_ptr(i) + j - row_start
            row_ptr(i+1) = row_start2
            row_end2 = row_ptr(i+2) - 1
            do k=row_end+1,row_end2
               k1 = row_start2 + k - (row_end+1)
               col_ind_0(k1) = col_ind_0(k)
               col_ind_0(k) = 0
            enddo
            goto 40
         endif
      enddo
40    continue
   enddo

   ! squeeze the last row
   i = neq
   row_start = row_ptr(i)
   row_end = row_ptr(i+1) - 1
   do j=row_start,row_end
      if(col_ind_0(j) .eq. 0) then
         row_start2 = row_ptr(i) + j - row_start
         row_ptr(i+1) = row_start2
         goto 50
      endif
   enddo
50 continue


   ! Find the longest row
   max_row_len = 0
   do i=1,neq
      row_start = row_ptr(i)
      row_end = row_ptr(i+1) - 1
      row_len = row_end - row_start + 1
      if (row_len .gt. max_row_len) max_row_len = row_len
   enddo

   if (debug .eq. 1) then
      write(ui,*) "csr_length: max_row_len = ", max_row_len
   endif

   ! Now we know nonz
   call integer_alloc_1d(col_ind, nonz, 'col_ind', errco, emsg); RETONERROR(errco)

   ! weird reverse labelleling because of reverse convention in this file
   col_ind(1:nonz) = col_ind_0(1:nonz)

   deallocate(col_ind_0)

   return
end
