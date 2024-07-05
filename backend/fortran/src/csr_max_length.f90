
subroutine csr_max_length (nel, n_ddl, neq, table_N_E_F, ineq, col_ptr, nonz)

   use numbatmod

   integer(8) nel, neq, n_ddl, nonz
   integer(8) table_N_E_F (14,nel)
   integer(8) ineq(3,n_ddl), col_ptr(neq+1)


   integer(8) i, k, iel, ind_ip, ip
   integer(8) k_copy1, k_copy2


   col_ptr = 0


   !  Determination of the bandwidths

   do  iel=1,nel
      do i=1,nddl_0_em
         ip = table_N_E_F(i,iel)
         do k=1,3
            ind_ip = ineq(k,ip)
            if (ind_ip .ne. 0) col_ptr(ind_ip) = col_ptr(ind_ip)+1
         enddo
      enddo
   enddo


   !  Compressed Row Storage (CRS): determine the row pointer

   k_copy1 = col_ptr(1)
   col_ptr(1) = 1
   do i=2,neq+1
      k_copy2 = col_ptr(i)
      col_ptr(i) = col_ptr(i-1) + 3*nddl_0_em + 3*(nddl_0_em-1)*(k_copy1-1)
      k_copy1 = k_copy2
   enddo

   nonz = col_ptr(neq+1) - 1

   return
end
