
subroutine csr_length (nel, n_ddl, neq,  &
   table_N_E_F, ineq, col_ind, row_ptr, nonz_max, &
   nonz, max_row_len, ipointer, int_max, debug)

   use numbatmod

   integer(8) nel, neq, n_ddl,  nonz_max, nonz
   integer(8) ipointer, int_max
   integer(8) table_N_E_F(14,nel)
   integer(8) ineq(3,n_ddl)
   integer(8) col_ind(*)
   integer(8) row_ptr(neq+1)
   integer(8) max_row_len


   integer alloc_stat
   integer(8), dimension(:), allocatable :: col_ind_0

   integer(8) i, j, k, k1, i_ddl, j_ddl
   integer(8) iel, ind_ip, ip, ind_jp, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2
   integer(8) ui, debug

   ui = stdout

   alloc_stat = 0

   allocate(col_ind_0(nonz_max), STAT=alloc_stat)
   if (alloc_stat /= 0) then
      write(*,*)
      write(*,*) "csr_length: ", "The allocation is unsuccessful"
      write(*,*) "alloc_stat = ", alloc_stat
      write(*,*) "Not enough memory for the array col_ind_0"
      write(*,*) "nonz_max = ", nonz_max
      write(*,*) "csr_length: Aborting..."
      stop
   endif



   !  do i=1,nonz_max
   !  col_ind_0(i) = 0
   !  enddo
   col_ind_0 = 0

!  Determination of the column indices

   nonz = 0
   do iel=1,nel
      do i=1,nddl_0_em
         ip = table_N_E_F(i,iel)
         do i_ddl=1,3
            ind_ip = ineq(i_ddl,ip)
            if (ind_ip .ne. 0) then
               row_start = row_ptr(ind_ip)
               row_end = row_ptr(ind_ip+1) - 1
               do j=1,nddl_0_em
                  jp = table_N_E_F(j,iel)
                  do j_ddl=1,3
                     ind_jp = ineq(j_ddl,jp)
                     if (ind_jp .ne. 0) then
!  Search if the entry (ind_ip,ind_jp) is already stored
                        do k=row_start,row_end
                           if(col_ind_0(k) .eq. 0) goto 20
                           if(col_ind_0(k) .eq. ind_jp) goto 30
                        enddo

                        print*, "csr_length: There is a problem!",&
                        &" Aborting..."
                        stop

20                      continue
!  No entry exists for (ind_ip,ind_jp); create new one
                        nonz = nonz + 1

                        if (nonz .gt. nonz_max) then
                           print*, "csr_length: nonz > nonz_max: ",&
                           &nonz .gt. nonz_max
                           print*, "csr_length: Aborting..."
                           stop
                        endif

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

   if (nonz .lt. nonz_max) then
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
40       continue
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
50    continue
   endif

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

   if ((ipointer+nonz) .gt. int_max) then
      write(ui,*) "csr_length: (ipointer+nonz) > int_max : ",&
      &(ipointer+nonz), int_max
      write(ui,*) "csr_length: nonz_max = ", nonz_max
      write(ui,*) "csr_length: increase the size of int_max"
      write(ui,*) "csr_length: Aborting..."
      stop
   else


!  Copy the local array col_ind_0 into col_ind
      !  do i=1,nonz
      !  col_ind(i) = col_ind_0(i)
      !  enddo
      col_ind(1:nonz) = col_ind_0(1:nonz)
   endif

   return
end
