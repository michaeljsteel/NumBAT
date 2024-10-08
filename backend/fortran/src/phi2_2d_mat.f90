

!  phi2_2d_mat evaluates a quadratic basis function (P2) and its derivative.
!  P2 basis function over the unit Triangle
!  phi(i) is the basis element i evaluated at vx=(x,y)
!  mat_grad(i,j) is Jacobian: x and y derivatives indexed by i, basis element by j
!    See lag_interpolation.nb
!


subroutine phi2_2d_mat (vx, phi, mat_grad)

   implicit none
   double precision vx(2), phi(6), mat_grad(2,6)
   double precision x0, y0
   integer(8) inode

   x0 = x(1)
   y0 = x(2)

   inode = 1
   phi(inode) = (-1 + x0 + y0)*(-1 + 2*x0 + 2*y0)
   mat_grad(1,inode) = -3 + 4*x0 + 4*y0
   mat_grad(2,inode) = -3 + 4*x0 + 4*y0

   inode = 2
   phi(inode) = x0*(-1 + 2*x0)
   mat_grad(1,inode) = -1 + 4*x0
   mat_grad(2,inode) = 0

   inode = 3
   phi(inode) = y0*(-1 + 2*y0)
   mat_grad(1,inode) = 0
   mat_grad(2,inode) = -1 + 4*y0

   inode = 4
   phi(inode) = -4*x0*(-1 + x0 + y0)
   mat_grad(1,inode) = -4*(-1 + 2*x0 + y0)
   mat_grad(2,inode) = -4*x0

   inode = 5
   phi(inode) = 4*x0*y0
   mat_grad(1,inode) = 4*y0
   mat_grad(2,inode) = 4*x0

   inode = 6
   phi(inode) = -4*y0*(-1 + x0 + y0)
   mat_grad(1,inode) = -4*y0
   mat_grad(2,inode) = -4*(-1 + x0 + 2*y0)



   return
end
