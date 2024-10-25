
!
!     P2 Lagrange Interpolation nodes for the unit triangle
!
!     unit triangle =  triangle whose vertices are:
!                         (0,0,0), (1,0,0), (0,1,0).
!

subroutine get_P2_node_locations (nnodes, xn)
    use numbatmod

   integer(8) nnodes
   double precision xn(2,P2_NODES_PER_EL)
   integer(8) i


   i = 1
   xn(1,i) = 0
   xn(2,i) = 0

   i = 2
   xn(1,i) = 1
   xn(2,i) = 0

   i = 3
   xn(1,i) = 0
   xn(2,i) = 1

   i = 4
   xn(1,i) = 0.5
   xn(2,i) = 0

   i = 5
   xn(1,i) = 0.5
   xn(2,i) = 0.5

   i = 6
   xn(1,i) = 0
   xn(2,i) = 0.5

end
