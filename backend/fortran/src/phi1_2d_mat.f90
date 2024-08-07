!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  phi1_2d_mat evaluates a linear basis function (P3) and its derivative.

!  P3 basis function over the unit Triangle

!  phi1_i = L^1(x) L^1(y)


subroutine phi1_2d_mat (x, phi, mat_grad)

   implicit none
   double precision x(2), phi(3), mat_grad(2,3)
   double precision x0, y0
   integer(8) inode

   x0 = x(1)
   y0 = x(2)

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
