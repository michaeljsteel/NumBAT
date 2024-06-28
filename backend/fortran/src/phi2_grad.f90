
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Note that if grad_i_0 is the gradient on the reference triangle,
!  then the gradient on the actual triangle is:
!  grad_i  = Transpose(mat_T)*grad_i0



subroutine phi2_grad(inode, nnodes, mat_jac, vec_grad)


   implicit none
   integer(8) inode, nnodes
   double precision vec_grad(2,nnodes)
   double precision mat_jac(2,2)

!  Local variables
   double precision c(2,2)
   integer(8) nnodes_0
   parameter (nnodes_0 = 6)
   double precision xel_0(2,nnodes_0)
   double precision phi_xi, phi_yi
   double precision phi0_xi, phi0_yi
   double precision x, y
   integer(8) i

!  integer(8) ui, debug

!  Coordinates (x,y)= xel_0(1..2,inode) of the P2 Lagrange interpolaion nodes
   xel_0(1,1) = 0
   xel_0(2,1) = 0

   xel_0(1,2) = 1
   xel_0(2,2) = 0

   xel_0(1,3) = 0
   xel_0(2,3) = 1

   xel_0(1,4) = 0.5d0
   xel_0(2,4) = 0

   xel_0(1,5) = 0.5d0
   xel_0(2,5) = 0.5d0

   xel_0(1,6) = 0
   xel_0(2,6) = 0.5d0

!  C = Tanspose[mat_jac]
   c(1,1) = mat_jac(1,1)
   c(2,2) = mat_jac(2,2)
   c(1,2) = mat_jac(2,1)
   c(2,1) = mat_jac(1,2)

   x = xel_0(1,inode)
   y = xel_0(2,inode)

   i = 1
!  x-derivative over the reference triangle
   phi0_xi = 4.0d0*(x+y) - 3.0d0

!  y-derivative over the reference triangle
   phi0_yi = 4.0d0*(x+y) - 3.0d0

!  x-derivative over the current triangle

   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi

!  y-derivative over the current triangle
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 2
   phi0_xi = 4.0d0*x - 1.0d0
   phi0_yi = 0.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 3
   phi0_xi = 0.0d0
   phi0_yi = 4.0d0*y - 1.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 4
   phi0_xi = 4.0d0*(1.0d0 - 2.0d0*x - y)
   phi0_yi = -4.0d0*x
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 5
   phi0_xi = 4.0d0*y
   phi0_yi = 4.0d0*x
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 6
   phi0_xi = -4.0d0*y
   phi0_yi = 4.0d0*(1.0d0 - x - 2.0d0*y)
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   return
end
