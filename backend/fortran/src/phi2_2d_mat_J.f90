
!  Evaluates the quadratic P2 polynomials on the six nodes a_i
!  This is the same as phi1_2d_mat.f90 without the derivatives
!  Despite the name, it is NOT a Jacobian

subroutine phi2_2d_mat_J (vx, phi)

   implicit none
   double precision vx(2), phi(6)
   double precision x0, y0
   integer(8) inode


   x0 = vx(1)
   y0 = vx(2)

   inode = 1
   phi(inode) = (-1 + x0 + y0)*(-1 + 2*x0 + 2*y0)

   inode = 2
   phi(inode) = x0*(-1 + 2*x0)

   inode = 3
   phi(inode) = y0*(-1 + 2*y0)

   inode = 4
   phi(inode) = -4*x0*(-1 + x0 + y0)

   inode = 5
   phi(inode) = 4*x0*y0

   inode = 6
   phi(inode) = -4*y0*(-1 + x0 + y0)


   return
end
