

!  Compute the gradient of the 10 P3 Lagrange polynomial basis function at one of the 6 P2 Lagrange nodes

!  Note that if grad_i_0 is the gradient on the reference triangle,
!  then the gradient on the actual triangle is:
!  grad_i  = Transpose(mat_T)*grad_i0

!  Finds the gradient of each basis function in the parametrised element
!  at node inode in [1..6]
!  Transformation using mat_jac, returned in vec_grad

subroutine phi3_grad_p2(inode, nnodes_P3, mat_jac, vec_grad)


   implicit none
   integer(8) inode, nnodes_P3
   double precision vec_grad(2,nnodes_P3)
   double precision mat_jac(2,2)

!  Local variables
   double precision c(2,2)
   integer(8), parameter :: nnodes_P2_0 = 6
   !double precision xel_0(2,nnodes_P2_0)

   double precision, dimension(2,nnodes_P2_0) :: xel_0 = &
      reshape( (/ 0.d0, 1.d0, 0.d0, 0.5d0, 0.5d0, 0.d0, &
      0.d0, 0.d0, 1.d0, 0.0d0, 0.5d0, 0.5d0 /), &
      shape(xel_0), order=(/2,1/))

   double precision phi_xi, phi_yi
   double precision phi0_xi, phi0_yi
   double precision x, y
   integer(8) i


!  Coordinates (x,y)= xel_0(1..2,inode) of the P2 Lagrange interpolaion nodes
   ! xel_0(1,1) = 0
   ! xel_0(2,1) = 0

   ! xel_0(1,2) = 1
   ! xel_0(2,2) = 0

   ! xel_0(1,3) = 0
   ! xel_0(2,3) = 1

   ! xel_0(1,4) = 0.5d0
   ! xel_0(2,4) = 0

   ! xel_0(1,5) = 0.5d0
   ! xel_0(2,5) = 0.5d0

   ! xel_0(1,6) = 0
   ! xel_0(2,6) = 0.5d0

!  C = Tanspose[mat_jac]
   c(1,1) = mat_jac(1,1)
   c(2,2) = mat_jac(2,2)
   c(1,2) = mat_jac(2,1)
   c(2,1) = mat_jac(1,2)

   x = xel_0(1,inode)
   y = xel_0(2,inode)

   i = 1
   !  x-derivative over the reference triangle
   phi0_xi = (-11 - 27*x**2 + x*(36 - 54*y) + 36*y - 27*y**2)/2.0d0

   !  y-derivative over the reference triangle
   phi0_yi = (-11 - 27*x**2 + x*(36 - 54*y) + 36*y - 27*y**2)/2.0d0

   !  x-derivative over the current triangle
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi

   !  y-derivative over the current triangle
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 2
   phi0_xi = 1 - 9*x + (27*x**2)/2.0d0
   phi0_yi = 0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 3
   phi0_xi = 0
   phi0_yi = 1 - 9*y + (27*y**2)/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 4
   phi0_xi = (9*(2 + 9*x**2 - 5*y + 3*y**2 + 2*x*(-5 + 6*y)))/2.0d0
   phi0_yi = (9*x*(-5 + 6*x + 6*y))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 5
   phi0_xi = (-9*(1 + 9*x**2 - y + x*(-8 + 6*y)))/2.0d0
   phi0_yi = (-9*x*(-1 + 3*x))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 6
   phi0_xi = (9*(-1 + 6*x)*y)/2.0d0
   phi0_yi = (9*x*(-1 + 3*x))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 7
   phi0_xi = (9*y*(-1 + 3*y))/2.0d0
   phi0_yi = (9*x*(-1 + 6*y))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 8
   phi0_xi = (-9*y*(-1 + 3*y))/2.0d0
   phi0_yi = (-9*(1 - 8*y + 9*y**2 + x*(-1 + 6*y)))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 9
   phi0_xi = (9*y*(-5 + 6*x + 6*y))/2.0d0
   phi0_yi = (9*(2 + 3*x**2 - 10*y + 9*y**2 + x*(-5 + 12*y)))/2.0d0
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi

   i = 10
   phi0_xi = -27*y*(-1 + 2*x + y)
   phi0_yi = -27*x*(-1 + x + 2*y)
   phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
   phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
   vec_grad(1,i) = phi_xi
   vec_grad(2,i) = phi_yi


   return
end
