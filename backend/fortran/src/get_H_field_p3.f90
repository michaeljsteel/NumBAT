
!  For a waveguide mode: compute the H-field on P2 nodes from the E-field using P2 and P3 data
!  Incoming Ez components are defined in NumBAT form: \hatEz=-i \beta Ez

! Applying the Maxwell's equations to the E-field of a waveguide mode with E=E0 exp(i (beta z-\omega t)), we get:
!    H= 1/(i mu0 omega) (zhat \cross dEt/dz - zhat \cross \gradt Ez + \gradt \cross Et)
!    = 1/(i mu0 omega) (zhat \cross (i\beta (Ex,Ey,0)) - zhat \cross (d Ez/dx, dEz/dy, 0)  + \gradt \cross Et)
!    H_x = [-i beta*E_y + D(E_z,y)] * gamma
!    H_y = [ i beta*E_x - D(E_z,x)] * gamma
!    H_z = [ D(E_y,x) - D(E_x,y)] * gamma
!  with gamma = 1/(i mu0 omega)


subroutine get_H_field_p3 (k_0, beta, mat_T, E_xyz_p2, E_z_p3, H_xyz_p2)

   use numbatmod

   double precision k_0, mat_T(2,2)
   complex(8) beta
   complex(8) E_xyz_p2(3,P2_NODES_PER_EL)
   complex(8) E_z_p3(10)
   complex(8) H_xyz_p2(3,P2_NODES_PER_EL)

   ! Locals

   double precision vec_grad_P2(2,P2_NODES_PER_EL)
   double precision vec_grad_P3(2,P3_NODES_PER_EL)
   double precision omega
   integer(8) nd_i, nd_j
   complex(8)  z_tmp1, z_tmp2

   H_xyz_p2 = C_ZERO

   ! i beta(-Ey, Ex, 0)
   do nd_i=1,P2_NODES_PER_EL
      H_xyz_p2(1, nd_i) = H_xyz_p2(1,nd_i) - C_IM_ONE * beta * E_xyz_p2(2,nd_i)
      H_xyz_p2(2, nd_i) = H_xyz_p2(2,nd_i) + C_IM_ONE * beta * E_xyz_p2(1,nd_i)
   enddo

   do nd_i=1,P2_NODES_PER_EL
      !       vec_grad_p2: contains the gradients of all 6 basis polynomials at the node nd_i
      call phi2_grad(nd_i, P2_NODES_PER_EL, mat_T, vec_grad_P2)
      call phi3_grad_p2(nd_i, P3_NODES_PER_EL, mat_T, vec_grad_P3)

      ! (D(E_z,y), -D(E_z,x), 0)
      do nd_j=1,P3_NODES_PER_EL
         H_xyz_p2(1,nd_i) = H_xyz_p2(1,nd_i) + vec_grad_P3(2,nd_j) * E_z_p3(nd_j)
         H_xyz_p2(2,nd_i) = H_xyz_p2(2,nd_i) -vec_grad_P3(1,nd_j) * E_z_p3(nd_j)
      enddo

      do nd_j=1,P2_NODES_PER_EL
         z_tmp1 = -vec_grad_P2(2,nd_j) * E_xyz_p2(1,nd_j)   ! - D(E_x,y)
         z_tmp2 =  vec_grad_P2(1,nd_j) * E_xyz_p2(2,nd_j)   !   D(E_y,x)
         H_xyz_p2(3,nd_i) = H_xyz_p2(3,nd_i) + z_tmp1+z_tmp2
      enddo

   enddo

   omega = k_0 * SI_C_SPEED

   H_xyz_p2 = H_xyz_p2 / (C_IM_ONE* omega * SI_MU_0)

end

