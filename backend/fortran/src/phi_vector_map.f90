
!  P2 basis function over the unit triangle

!  Quadratic basis function = P2 * Grad P1

!  vector_elt_map(1,j,i) = k : number on data to be stored: if k=3 only one gradient will be used; k=4 => 2 gradients
!  vector_elt_map(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  vector_elt_map(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial phi_n
!  vector_elt_map(4,j,i)     : it will be used only if k=4

!  The basis functions are stored in the 3rd index i as : centre, edges 1(4), 2(5), 3(6)
!                                 in the 2nd index j as 3 functions on each edge (for i=1), and along each edge (for i=2..4)

subroutine make_vector_elt_map (nod_el, vector_elt_map)

   use numbatmod

   integer(8) nod_el(P2_NODES_PER_EL), vector_elt_map(4,3,N_ETY_TRANSVERSE)

   integer ety, bf, dof, j
   integer(8) list_end(2,3), j2
   integer(8) ls_n(3) !, ls_n_sorted(3)


   integer(8) edge_end, elo, ehi, iedge


!  Endpoints of the 6 edges (mid-point) of the reference triangle

   iedge = 1
   list_end(1,iedge) = 1
   list_end(2,iedge) = 2

   iedge = 2
   list_end(1,iedge) = 2
   list_end(2,iedge) = 3

   iedge = 3
   list_end(1,iedge) = 1
   list_end(2,iedge) = 3


   !  scan the element face
   ety=1

   !  The mid-edge nodes of the face
   !  scan the mid-edge nodes of the face

   do bf=1,3

      vector_elt_map(1,bf,ety) = 3        !  number on data to be stored

      vector_elt_map(2,bf,ety) = bf+3      !  the mid-edge number

      !  give the node opposite to the mid-edge node (j+3)
      j2 = modulo(bf+2,3)
      if( j2 .eq. 0 ) j2 = 3

      vector_elt_map(3,bf,ety) = j2
      ! second gradient not needed for this function
      vector_elt_map(4,bf,ety) = 0
   enddo

   !  scan the 3 element edges
   do ety=2,4
      !  2 end-point basis vectors are attached to the edge i
      !  scan the end nodes of the edge

      iedge = ety-1

      ! Find the indices corresponding to the nodes at the ends of the current edge
      ! They need to be sorted by the absolute value of the nodes in the overal mesh table (via nod_el)
      do j=1,2
         edge_end = list_end(j,iedge)
         ls_n(j) = nod_el(edge_end)
      enddo

      elo=1
      ehi=2
      ! swap if absolute node ordering is reversed
      if (ls_n(1) .gt. ls_n(2)) then
         elo=2
         ehi=1
      endif

      elo = list_end(elo,iedge)
      ehi = list_end(ehi,iedge)

      bf = 1
      vector_elt_map(1,bf,ety) = 3
      vector_elt_map(2,bf,ety) = elo
      vector_elt_map(3,bf,ety) = ehi
      vector_elt_map(4,bf,ety) = 0

      bf = 2
      vector_elt_map(1,bf,ety) = 3
      vector_elt_map(2,bf,ety) = ehi
      vector_elt_map(3,bf,ety) = elo
      vector_elt_map(4,bf,ety) = 0

      bf = 3
      vector_elt_map(1,bf,ety) = 4
      vector_elt_map(2,bf,ety) = iedge+3   !  add 3 to get the correct edge number
      vector_elt_map(3,bf,ety) = elo
      vector_elt_map(4,bf,ety) = ehi

   enddo

end
