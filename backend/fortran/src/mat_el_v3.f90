subroutine mat_el_v3 (basfuncs, el_nds_xy, beta, c_tensor_el, rho_el, mat_K, mat_M, nberr)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   type(NBError) nberr

   double precision el_nds_xy(2,6)
   complex(8) beta
   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el

   !  Local variables
   type(BasisFunctions) basfuncs


   ! double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
   ! double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
   ! double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
   complex(8) z_tmp1, z_tmp2
   complex(8) z_mat_xyz(6,6,3,3), z_beta, z_tensor, t_zmat_xyz
   integer(8) i, j, i_p, j_p, i_xyz,  j_xyz
   integer(8) i_u_xyz,  j_u_xyz, i_ind, j_ind
   integer(8) S_index(3,3)

   !  Compute the Affine mappings from the current triangle to the
   !  reference unit triangle.
   !  Integration will be performed on the reference unit triangle


   !  S_index(i_xyz,i_u_xyz): index of S_{i_xyz,i_u_xyz} in the Voigt notation
   !  i_xyz represents the x,y, or z-derivative
   !  i_u_xyz represents the x,y, or z-field component

   S_index(1,1) = 1 !  S_xx => 1
   S_index(2,1) = 6 !  S_yx => 6
   S_index(3,1) = 5 !  S_zx => 5

   S_index(1,2) = 6 !  S_xy => 6
   S_index(2,2) = 2 !  S_yy => 2
   S_index(3,2) = 4 !  S_zy => 4

   S_index(1,3) = 5 !  S_xz => 5
   S_index(2,3) = 4 !  S_yz => 4
   S_index(3,3) = 3 !  S_zz => 3



   !  z_mat_xyz: contains the overlap integrals of the x,y and z-derivatives
   z_mat_xyz = C_ZERO
   do i=1,P2_NODES_PER_EL
      do j=1,P2_NODES_PER_EL

         do i_xyz=1,3
            do j_xyz=1,3

               if (i_xyz == 1 .and. j_xyz == 1) then             ! u_i_x^*  u_j_x
                  t_zmat_xyz =  basfuncs%p2x_p2x(i,j)

               elseif (i_xyz == 1 .and. j_xyz == 2) then         ! u_i_x^*  u_j_y
                  t_zmat_xyz = basfuncs%p2x_p2y(i,j)

               elseif (i_xyz == 1 .and. j_xyz == 3) then         ! u_i_x^*  u_j_z
                  t_zmat_xyz = C_IM_ONE* beta * basfuncs%p2_p2x(j,i)

               elseif (i_xyz == 2 .and. j_xyz == 1) then         ! u_i_y^* u_j_x
                  t_zmat_xyz = basfuncs%p2x_p2y(j,i)

               elseif (i_xyz == 2 .and. j_xyz == 2) then         ! u_i_y^* u_j_y
                  t_zmat_xyz = basfuncs%p2y_p2y(i,j)

               elseif (i_xyz == 2 .and. j_xyz == 3) then         ! u_i_y^* u_j_z
                  t_zmat_xyz = C_IM_ONE* beta* basfuncs%p2_p2y(j,i)

               elseif (i_xyz == 3 .and. j_xyz == 1) then         ! u_i_z^* u_j_x
                  t_zmat_xyz = - C_IM_ONE* beta * basfuncs%p2_p2x(i,j)

               elseif (i_xyz == 3 .and. j_xyz == 2) then         ! u_i_z^* u_j_y
                  t_zmat_xyz = - C_IM_ONE* beta * basfuncs%p2_p2y(i,j)

               elseif (i_xyz == 3 .and. j_xyz == 3) then         ! u_i_z^* u_j_z
                  t_zmat_xyz = beta**2 * basfuncs%p2_p2(i,j)

               endif

               z_mat_xyz(i,j,i_xyz,j_xyz) = t_zmat_xyz

            enddo
         enddo

      enddo
   enddo

   mat_K = C_ZERO
   mat_M = C_ZERO

   !=================  Construction of the matrix mat_M =================

   !  Integral [rho * P(i) * P(i)]
   do i=1,P2_NODES_PER_EL

      do i_xyz=1,3  !  The components x, y and z
         i_p = 3*(i-1) + i_xyz

         do j=1,P2_NODES_PER_EL
            j_xyz = i_xyz
            j_p = 3*(j-1) + j_xyz
            mat_M(i_p,j_p) = mat_M(i_p,j_p) + basfuncs%p2_p2(i,j) * rho_el
         enddo

      enddo
   enddo

   !=================  Construction of the matrix mat_K =================
   !  Integral [K_{ij} = Gradient_s(conjg(P_k(i))) x c_tensor x Gradient_s(P_k(j))], where k=x,y,z
   !  Reference: see Eqs. (7) and (8) in:
   !  A.-C. Hladky-Hennion
   !  "Finite element analysis of the propagation of acoustic waves in waveguides,"
   !  Journal of Sound and Vibration, vol. 194, no. 2, pp. 119-136, 1996.

   do i=1,P2_NODES_PER_EL
      !  Components of the displacement vector

      do i_u_xyz=1,3
         i_p = 3*(i-1) + i_u_xyz

         do j=1,P2_NODES_PER_EL
            !  Components of the displacement vector
            do j_u_xyz=1,3
               j_p = 3*(j-1) + j_u_xyz
               !  Derivatives

               do i_xyz=1,3
                  i_ind = S_index(i_xyz,i_u_xyz)

                  !  Derivatives
                  do j_xyz=1,3
                     j_ind = S_index(j_xyz,j_u_xyz)
                     z_tensor = c_tensor_el(i_ind,j_ind)
                     z_tmp1 = z_mat_xyz(i,j,i_xyz,j_xyz)

                     if (i_u_xyz == 3) then
                        z_tmp1 = -C_IM_ONE* z_tmp1
                     endif

                     if (j_u_xyz == 3) then
                        z_tmp1 = C_IM_ONE* z_tmp1
                     endif

                     mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1 * z_tensor

                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo


end

