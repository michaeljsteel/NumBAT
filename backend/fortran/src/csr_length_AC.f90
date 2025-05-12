
#include "numbat_decl.h"

subroutine csr_length_AC (mesh, cscmat, &
   n_nonz_max,&
  n_nonz, max_row_len, nberr)

   use numbatmod
   use alloc

   use class_Mesh
   use class_SparseCSC_AC

   type(MeshAC) mesh
   type(SparseCSC_AC) cscmat
   type(NBError) nberr

   integer(8) n_nonz_max,n_nonz




   integer(8) max_row_len


   integer(8), dimension(:), allocatable :: row_ind_tmp

   integer(8), parameter :: N_ENTITY_PER_EL_AC = 6

   integer(8) i, j, i_nd, j_nd, k, k1, i_dof, j_dof
   integer(8) iel, ind_ip, ip, ind_jp, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2
   integer(8) ui



   ui = stdout

   call integer_alloc_1d(row_ind_tmp, n_nonz_max, 'row_ind_tmp', nberr); RET_ON_NBERR(nberr)

   row_ind_tmp = 0

   ! This code was converted from one in CSR format
   !  Determination of the row indices

  n_nonz = 0
   do iel=1,mesh%n_msh_elts

      do i_nd=1,N_ENTITY_PER_EL_AC
         ip = mesh%m_elnd_to_mshpt(i_nd,iel)

         do i_dof=1,3
            ind_ip = cscmat%m_eqs(i_dof,ip)
            if (ind_ip .ne. 0) then
               row_start = cscmat%v_col_ptr(ind_ip)
               row_end = cscmat%v_col_ptr(ind_ip+1) - 1

               do j_nd=1,N_ENTITY_PER_EL_AC
                  jp = mesh%m_elnd_to_mshpt(j_nd,iel)

                  do j_dof=1,3
                     ind_jp = cscmat%m_eqs(j_dof,jp)
                     if (ind_jp .ne. 0) then
                        ! Search if the entry (ind_ip,ind_jp) is already stored
                        do k=row_start,row_end
                           if(row_ind_tmp(k) .eq. 0) goto 20
                           if(row_ind_tmp(k) .eq. ind_jp) goto 30
                        enddo

                        print*, "csr_length_AC: There is a problem!", " Aborting..."
                        stop

20                      continue

                        !  No entry exists for (ind_ip,ind_jp); create new one
                       n_nonz =n_nonz + 1
                        if (n_nonz .gt. n_nonz_max) then
                           print*, "csr_length_AC:n_nonz > n_nonz_max: ",&
                           n_nonz .gt. n_nonz_max
                           print*, "csr_length_AC: Aborting..."
                           stop
                        endif

                        row_ind_tmp(k) = ind_jp
30                      continue

                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo


   ! squeeze away the zero entries
   ! added so as to handle more type of domains/meshes

   if (n_nonz .lt. n_nonz_max) then
      do i=1,cscmat%n_dof-1
         row_start = cscmat%v_col_ptr(i)
         row_end = cscmat%v_col_ptr(i+1) - 1
         do j=row_start,row_end
            if(row_ind_tmp(j) .eq. 0) then
               row_start2 = cscmat%v_col_ptr(i) + j - row_start
               cscmat%v_col_ptr(i+1) = row_start2
               row_end2 = cscmat%v_col_ptr(i+2) - 1
               do k=row_end+1,row_end2
                  k1 = row_start2 + k - (row_end+1)
                  row_ind_tmp(k1) = row_ind_tmp(k)
                  row_ind_tmp(k) = 0
               enddo
               goto 40
            endif
         enddo
40       continue
      enddo

      i = cscmat%n_dof
      row_start = cscmat%v_col_ptr(i)
      row_end = cscmat%v_col_ptr(i+1) - 1
      do j=row_start,row_end
         if(row_ind_tmp(j) .eq. 0) then
            row_start2 = cscmat%v_col_ptr(i) + j - row_start
            cscmat%v_col_ptr(i+1) = row_start2
            goto 50
         endif
      enddo
50    continue
   endif


   return
   max_row_len = 0
   do i=1,cscmat%n_dof
      row_start = cscmat%v_col_ptr(i)
      row_end = cscmat%v_col_ptr(i+1) - 1
      row_len = row_end - row_start + 1
      if (row_len .gt. max_row_len) max_row_len = row_len
   enddo

   write(*,*) 'csclength', max_row_len




   ! Now we know n_nonz
   !call integer_alloc_1d(cscmat%v_row_ind, n_nonz, 'this%v_row_ind', nberr); RET_ON_NBERR(nberr)

   !cscmat%v_row_ind(1:n_nonz) = row_ind_tmp(1:n_nonz)



end
