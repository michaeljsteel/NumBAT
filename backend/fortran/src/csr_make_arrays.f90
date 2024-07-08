#include "numbat_decl.h"

subroutine make_csr_arrays(n_msh_el, n_ddl, neq, table_N_E_F, &
   m_eqs, nonz, v_row_ind, v_col_ptr, debug, errco, emsg)

   use numbatmod
   use alloc
   use nbinterfacesb

   integer(8) n_msh_el, n_ddl, neq, nonz
   integer(8) table_N_E_F(14,n_msh_el)
   integer(8) m_eqs(3,n_ddl)

   integer(8), dimension(:), allocatable, intent(inout) :: v_row_ind
   integer(8), dimension(:), allocatable, intent(inout) :: v_col_ptr

   integer(8) debug

   integer errco
   character(len=EMSG_LENGTH) emsg

   ! ------------------------------------------

   integer(8) nonz_max, max_row_len
   integer(8), dimension(:), allocatable  :: iwork

   ! ------------------------------------------



   call integer_alloc_1d(v_col_ptr, neq+1, 'v_col_ptr', errco, emsg); RETONERROR(errco)

   call csr_max_length (n_msh_el, n_ddl, neq, table_N_E_F, &
      m_eqs, v_col_ptr, nonz_max)


   ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
   ! length of v_row_ind is determined inside csr_length and so allocated there
   call csr_length (n_msh_el, n_ddl, neq,  table_N_E_F, m_eqs, &
      v_row_ind, v_col_ptr, nonz_max, nonz, max_row_len, debug, errco, emsg)
   RETONERROR(errco)

   call integer_alloc_1d(iwork, 3*n_ddl, 'iwork', errco, emsg);
   RETONERROR(errco)

   call sort_csr (neq, nonz, max_row_len, v_row_ind, v_col_ptr,  iwork)


end subroutine

! row/col names seem backward
! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine csr_length (n_msh_el, n_ddl, neq,  &
   table_N_E_F, m_eqs, &
   col_ind, row_ptr, &  ! these names are swtiched from the call, but matched to the weird reverse naming in this file
   nonz_max, nonz, max_row_len, debug, errco, emsg)

   use numbatmod
   use alloc

   integer(8) n_msh_el, n_ddl, neq
   integer(8) table_N_E_F(14,n_msh_el)
   integer(8) m_eqs(3,n_ddl)


   integer(8), dimension(:), allocatable, intent(inout) :: col_ind
   integer(8), dimension(:) :: row_ptr(neq+1)

   integer(8) nonz_max, nonz, max_row_len

   integer(8) debug

   integer errco
   character(len=EMSG_LENGTH) emsg

   ! --------------------------------------------

   integer(8), dimension(:), allocatable :: col_ind_0

   integer(8) i, j, k, k1, i_ddl, j_ddl
   integer(8) iel, ind_ip, ip, ind_jp, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2, ui


   ui = stdout


   call integer_alloc_1d(col_ind_0, nonz_max, 'col_ind_0', errco, emsg); RETONERROR(errco)

   col_ind_0 = 0

   !  Determination of the column indices

   nonz = 0
   do iel=1,n_msh_el
      do i=1,nddl_0_em
         ip = table_N_E_F(i,iel)
         do i_ddl=1,3
            ind_ip = m_eqs(i_ddl,ip)

            if (ind_ip .ne. 0) then
               row_start = row_ptr(ind_ip)
               row_end = row_ptr(ind_ip+1) - 1
               do j=1,nddl_0_em
                  jp = table_N_E_F(j,iel)
                  do j_ddl=1,3
                     ind_jp = m_eqs(j_ddl,jp)

                     if (ind_jp .ne. 0) then
!  Search if the entry (ind_ip,ind_jp) is already stored
                        do k=row_start,row_end
                           if(col_ind_0(k) .eq. 0) goto 20
                           if(col_ind_0(k) .eq. ind_jp) goto 30
                        enddo

                        ! bail out
                        emsg = "csr_length: There is a problem with row/col indexing!"
                        errco = NBERROR_118
                        return

20                      continue
!  No entry exists for (ind_ip,ind_jp); create new one
                        nonz = nonz + 1

                        col_ind_0(k) = ind_jp

30                      continue
                     endif
                  enddo
               enddo
            endif
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

   ! weird rreverse labelleling because of reverse convention in this file
   col_ind(1:nonz) = col_ind_0(1:nonz)

   deallocate(col_ind_0)

   return
end
