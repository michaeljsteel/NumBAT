
subroutine csr_make_col_ptr_loose_AC (nel, npt, n_dof, nnodes, &
   elnd_to_mshpt, in_dof, lb, nonz)

   use numbatmod

   implicit none
   integer(8),  intent(in) :: nnodes
   integer(8), intent(in) :: nel, n_dof, npt
   integer(8), intent(in) ::  elnd_to_mshpt (nnodes,nel)
   integer(8), intent(in) :: in_dof(3,npt)
   integer(8), intent(out) :: lb(n_dof+1)
   integer(8), intent(out) :: nonz


   !  Local variables
   integer(8), parameter :: nddl_0 = 6 !  Different to em case. Different FEM?

   integer(8) :: i, k, iel, ind_ip, ip
   integer(8) :: k_copy1, k_copy2

   call assert_or_die(nnodes == 6,  &
      "csr_make_col_ptr_loose problem: nnodes = " // int_2_str(nnodes) // "." // &
      new_line('a') // " n_nodes should be equal to 6 !" // new_line('a') // &
      " Aborting...",  20_8)

   !do i=1,n_dof+1
   !lb(i) = 0
   !enddo
   lb = 0


   !  Determination of the bandwidths


   do iel=1,nel
      do i=1,nddl_0
         ip = elnd_to_mshpt(i,iel)
         do k=1,3
            ind_ip = in_dof(k,ip)
            if (ind_ip .ne. 0) lb(ind_ip) = lb(ind_ip)+1
         enddo
      enddo
   enddo

   nonz = 0
   do i=1,n_dof
      nonz = nonz + 3*nddl_0 + 3*(nddl_0-1)*(lb(i)-1)
   enddo


   !  Compressed Row Storage (CRS): determine the row pointer


   k_copy1 = lb(1)
   lb(1) = 1
   do i=2,n_dof+1
      k_copy2 = lb(i)
      lb(i) = lb(i-1) + 3*nddl_0 + 3*(nddl_0-1)*(k_copy1-1)
      k_copy1 = k_copy2
   enddo

   nonz = lb(n_dof+1) - 1


   return
end
