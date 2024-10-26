!  P2 vector basis functions over the reference unit triangle

!  Compute:
!  a quadratic vector basis function (vec_phi = P2 * Grad P1) and its transverse curl (curlt_phi)

! vector_elt_map indexes the pieces of the basis function
!  vector_elt_map(1,j,i) = k : number of gradients needed: if k=3 only one gradient will be used; k=4 => 2 gradients
!  vector_elt_map(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  vector_elt_map(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial eta_n
!  vector_elt_map(4,j,i)     : gradient of the second P1 polynomial eta_n2 if needed


! Generate the vector basis function for the bf_j bf_j of transverse entity ety_trans and its transverse curl (purely z-component)
! in terms of the P1 and P2 scalar functions

! formulas
! ety_trans=1, face element
!    3 edge functions at P2 nodes 4, 5, 6
!       j=1, node 4: phi_4 Grad eta_3
!       j=2, node 5: phi_5 Grad eta_1
!       j=3, node 6: phi_6 Grad eta_2

! ety_trans=2, edge 1, end nodes 1&2
!       j=1, : phi_1 Grad phi_2
!       j=2, : phi_2 Grad phi_1
!       j=3, : phi_4 (Grad phi_1-Grad phi_2)

! ety_trans=3, edge 2, end nodes 2&3
!       j=1, : phi_2 Grad phi_3
!       j=2, : phi_3 Grad phi_2
!       j=3, : phi_5 (Grad phi_2-Grad phi_3)

! ety_trans=4, edge 3, end nodes 3&1
!       j=1, : phi_3 Grad phi_1
!       j=2, : phi_1 Grad phi_3
!       j=3, : phi_6 (Grad phi_1-Grad phi_3)

subroutine evaluate_vector_elts(bf_j, ety_trans, vector_elt_map, phi_P2_ref,&
   gradt_P1_act, gradt_P2_act, vec_phi, curlt_phi)

   use numbatmod

   integer(8) bf_j, ety_trans
   integer(8) vector_elt_map(4,3,N_ETY_TRANSVERSE)
   double precision phi_P2_ref(P2_NODES_PER_EL)
   double precision gradt_P1_act(2,3), gradt_P2_act(2,P2_NODES_PER_EL)
   double precision vec_phi(2), curlt_phi

   integer debug
   !  Local variables
   integer(8) k, m, n1, n2
   double precision grad_p1(2), grad_p2(2), phi


   debug = 0

   if (debug .gt. 0) then
      if (k .eq. 4 .and. n2 .lt. 1) then
         write(*,*) "basis_vec: problem n2 < 1 for k = 4 "
         write(*,*) "basis_vec: n2 should >= 1 for k=4 !"
         write(*,*) "basis_vec: k, m, n1, n2 = ", k, m, n1, n2
         write(*,*) "basis_vec: Aborting..."
         stop
      endif
   endif

   k  = vector_elt_map(1, bf_j, ety_trans)
   m  = vector_elt_map(2, bf_j, ety_trans)
   n1 = vector_elt_map(3, bf_j, ety_trans)
   n2 = vector_elt_map(4, bf_j, ety_trans)



   if (k .eq. 3) then

      phi = phi_P2_ref(m)
      grad_p2 = gradt_P2_act(:,m)
      grad_p1 = gradt_P1_act(:,n1)
      vec_phi = phi * grad_p1


   elseif (k .eq. 4) then

      phi = phi_P2_ref(m)
      grad_p2 = gradt_P2_act(:,m)
      grad_p1 = gradt_P1_act(:,n1) - gradt_P1_act(:,n2)
      vec_phi = phi * grad_p1


   endif

   !  Curl_t E = Det( grad_p2,  grad_p1)
   curlt_phi = grad_p2(1)*grad_p1(2) - grad_p2(2)*grad_p1(1)

end
