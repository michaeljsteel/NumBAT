


!  P2 basis function over the unit Tetrahedron

!  Compute:
!  a quadradic basis function (vec_phi = P2 * Grad P1) and
!  and its transverse curl (curl_t_phi)

!  basis_list(1,j,i) = k : number on data to be stored: if k=3 only one gradient will be used; k=4 => 2 gradients
!  basis_list(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  basis_list(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial phi_n
!  basis_list(4,j,i)     : it will be used only if k=4


! generate the vector function i_ddl and its curl in terms of the P1 and P2 functions
subroutine basis_vec (i_eq, i_ddl, basis_list, p2_list,&
grad_p1_mat, grad_p2_mat, vec_phi, curl_t_phi)

   use numbatmod

   integer(8), parameter :: nnodes = 6
   !integer(8), parameter :: nddl_t=4
   integer(8), parameter :: dimm=2

   integer(8) i_eq, i_ddl



   integer(8) basis_list(4,3,nddl_t)
   double precision p2_list(nnodes)
   double precision grad_p1_mat(dimm,3), grad_p2_mat(dimm,nnodes)
   double precision vec_phi(dimm), curl_t_phi

   !  Local variables
   integer(8) i, k, m, n1, n2
   double precision grad_p1(dimm), grad_p2(dimm), phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   k  = basis_list(1, i_eq, i_ddl)
   m  = basis_list(2, i_eq, i_ddl)
   n1 = basis_list(3, i_eq, i_ddl)
   n2 = basis_list(4, i_eq, i_ddl)

   if (k .eq. 3) then

      phi = p2_list(m)
      do i=1,dimm
         grad_p2(i) = grad_p2_mat(i,m)
         grad_p1(i) = grad_p1_mat(i,n1)
         vec_phi(i) = phi * grad_p1(i)
      enddo

      elseif (k .eq. 4) then

         if (n2 .lt. 1) then
         write(*,*) "basis_vec: problem n2 < 1 for k = 4 "
         write(*,*) "basis_vec: n2 should >= 1 for k=4 !"
         write(*,*) "basis_vec: k, m, n1, n2 = ", k, m, n1, n2
         write(*,*) "basis_vec: Aborting..."
         stop
      endif

      phi = p2_list(m)
      do i=1,dimm
         grad_p2(i) = grad_p2_mat(i,m)
         grad_p1(i) = grad_p1_mat(i,n1) - grad_p1_mat(i,n2)
         vec_phi(i) = phi * grad_p1(i)
      enddo

   else
      write(*,*) "basis_vec: no action is defined when k = ", k
      write(*,*) "basis_vec: k should be equal to 3 or 4"
      write(*,*) "basis_vec: Aborting..."
      stop
   endif

   !  Curl_t E = Det( grad_p2,  grad_p1)
   curl_t_phi = grad_p2(1)*grad_p1(2) - grad_p2(2)*grad_p1(1)

   return
end
