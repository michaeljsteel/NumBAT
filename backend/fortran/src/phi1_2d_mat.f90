!  phi1_2d_mat evaluates a linear P1 basis function and its derivative
!  at a single point vx
!
!  phi1_i = L^1(x) L^1(y)
!  See lag_interpolation.nb

subroutine phi1_2d_mat (vx, phi, mat_grad)

   implicit none
   double precision vx(2), phi(3), mat_grad(2,3)
   double precision x0, y0
   integer(8) inode

   x0 = vx(1)
   y0 = vx(2)

   inode = 1
   phi(inode) = 1 - x0 - y0
   mat_grad(1,inode) = -1
   mat_grad(2,inode) = -1

   inode = 2
   phi(inode) = x0
   mat_grad(1,inode) = 1
   mat_grad(2,inode) = 0

   inode = 3
   phi(inode) = y0
   mat_grad(1,inode) = 0
   mat_grad(2,inode) = 1


   return
end
