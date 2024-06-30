
subroutine csr_max_length (nel, n_ddl, neq, table_N_E_F, ineq, lb, nonz)

   use numbatmod

   integer(8) nel, neq, n_ddl, nonz
   integer(8) table_N_E_F (14,nel)
   integer(8) ineq(3,n_ddl), lb(neq+1)


   integer(8) i, k, iel, ind_ip, ip
   integer(8) k_copy1, k_copy2


   lb = 0

   !  Determination of the bandwidths

   do  iel=1,nel
      do i=1,nddl_0_em
         ip = table_N_E_F(i,iel)
         do k=1,3
            ind_ip = ineq(k,ip)
            if (ind_ip .ne. 0) lb(ind_ip) = lb(ind_ip)+1
         enddo
      enddo
   enddo

   !  TODO: This block seems pointless. Overridden by later part
   nonz = 0
   do i=1,neq
      nonz = nonz + 3*nddl_0_em + 3*(nddl_0_em-1)*(lb(i)-1)
   enddo


   !  Compressed Row Storage (CRS): determine the row pointer

   k_copy1 = lb(1)
   lb(1) = 1
   do i=2,neq+1
      k_copy2 = lb(i)
      lb(i) = lb(i-1) + 3*nddl_0_em + 3*(nddl_0_em-1)*(k_copy1-1)
      k_copy1 = k_copy2
   enddo

   nonz = lb(neq+1) - 1


   return
end
