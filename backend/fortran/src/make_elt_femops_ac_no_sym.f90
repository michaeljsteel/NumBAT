! Construct the FEM matrix elements for the elastic problem at a particular element
! We are solving
! \mathrm{K} \vecu_h = \Omega^2 \mathrm{M} \vecu_h ,
! where the *stiffness* K and *mass* M matrices are defined as
!   K_{lm} = & \int_A (\nabla_s \vecg_l)^* : (\bar{c} : \nabla_s \vecg_m) \, \dA
!   M_{lm} = & \int_A \rho  \vecg_l^* \cdot  \vecg_m \, \dA
!
! basfuncs contains the basis functions and overlap integrals evaluated at the current element
! stiffness_C_IJ and rho_el are the stiffness and density for this element
! The K and M matrices are 18x18 from the 18 local degrees of freedom:
!   6 nodes each with 3 displacement components on each triangle

subroutine make_elt_femops_ac_no_sym (basfuncs, beta, stiffness_C_IJ, rho_el, mat_K, mat_M)

   use numbatmod
   use class_TriangleIntegrators
   use class_BasisFunctions

   complex(8) beta
   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) stiffness_C_IJ(6,6), rho_el

   !  Local variables
   type(BasisFunctions) basfuncs

   complex(8) t_G_ijbc, t_K_ij
   complex(8) Gmat_ij_bc(6,6,3,3), t_stiffC_IJ, tG_ij_bc
   integer(8) vnd_i, und_j, i_p, j_p, b_xyz,  c_xyz
   integer(8) v_xyz_i,  u_xyz_j, vgt_I, vgt_J
   integer(8) S_index(3,3)

   !  Compute the Affine mappings from the current triangle to the
   !  reference unit triangle.
   !  Integration will be performed on the reference unit triangle


   !  S_index(b_xyz,v_xyz_i): index of S_{b_xyz,v_xyz_i} in the Voigt notation
   !  b_xyz represents the x,y, or z-derivative
   !  v_xyz_i represents the x,y, or z-field component

   S_index(1,1) = 1 !  S_xx => 1
   S_index(2,1) = 6 !  S_yx => 6
   S_index(3,1) = 5 !  S_zx => 5

   S_index(1,2) = 6 !  S_xy => 6
   S_index(2,2) = 2 !  S_yy => 2
   S_index(3,2) = 4 !  S_zy => 4

   S_index(1,3) = 5 !  S_xz => 5
   S_index(2,3) = 4 !  S_yz => 4
   S_index(3,3) = 3 !  S_zz => 3


      mat_K = C_ZERO
      mat_M = C_ZERO

      !=================  Construction of the matrix mat_M =================
      ! See docs chap 9.
      ! M_{i,sig,j,tau} = delta_{sig,tau} \int_A \rho  g_i g_j \dx \dy
      !  Integral [rho * P(i) * P(i)]
      do vnd_i=1,P2_NODES_PER_EL

         do b_xyz=1,3  !  The components x, y and z
            i_p = 3*(vnd_i-1) + b_xyz

            c_xyz = b_xyz            ! only diagonal component entries are nonzero
            do und_j=1,P2_NODES_PER_EL
               j_p = 3*(und_j-1) + c_xyz
               mat_M(i_p,j_p) = mat_M(i_p,j_p) + basfuncs%p2_p2(vnd_i, und_j) * rho_el
            enddo

         enddo
      enddo


   ! Construction of the P2 derivative overlaps
   ! This is exactly the matrix G_ijbc in chap 9

   Gmat_ij_bc = C_ZERO
   do vnd_i=1,P2_NODES_PER_EL      ! g_i
      do und_j=1,P2_NODES_PER_EL   ! g_j

         do b_xyz=1,3              ! deriv of g_i
            do c_xyz=1,3           ! deriv of g_j

               if (b_xyz == 1 .and. c_xyz == 1) then             ! v_i_x^*  u_j_x
                  tG_ij_bc =  basfuncs%p2x_p2x(vnd_i, und_j)

               elseif (b_xyz == 1 .and. c_xyz == 2) then         ! v_i_x^*  u_j_y
                  tG_ij_bc = basfuncs%p2x_p2y(vnd_i,und_j)

               elseif (b_xyz == 1 .and. c_xyz == 3) then         ! v_i_x^*  u_j_z
                  tG_ij_bc = C_IM_ONE* beta * basfuncs%p2_p2x(und_j,vnd_i)

               elseif (b_xyz == 2 .and. c_xyz == 1) then         ! v_i_y^* u_j_x
                  tG_ij_bc = basfuncs%p2x_p2y(und_j,vnd_i)

               elseif (b_xyz == 2 .and. c_xyz == 2) then         ! v_i_y^* u_j_y
                  tG_ij_bc = basfuncs%p2y_p2y(vnd_i,und_j)

               elseif (b_xyz == 2 .and. c_xyz == 3) then         ! v_i_y^* u_j_z
                  tG_ij_bc = C_IM_ONE* beta* basfuncs%p2_p2y(und_j,vnd_i)

               elseif (b_xyz == 3 .and. c_xyz == 1) then         ! v_i_z^* u_j_x
                  tG_ij_bc = - C_IM_ONE* beta * basfuncs%p2_p2x(vnd_i,und_j)

               elseif (b_xyz == 3 .and. c_xyz == 2) then         ! v_i_z^* u_j_y
                  tG_ij_bc = - C_IM_ONE* beta * basfuncs%p2_p2y(vnd_i,und_j)

               elseif (b_xyz == 3 .and. c_xyz == 3) then         ! v_i_z^* u_j_z
                  tG_ij_bc = beta**2 * basfuncs%p2_p2(vnd_i,und_j)

               endif

               Gmat_ij_bc(vnd_i, und_j, b_xyz,c_xyz) = tG_ij_bc

            enddo
         enddo

      enddo
   enddo

   !=================  Construction of the matrix mat_K =================
   !  Integral [K_{ij} = Gradient_s(conjg(P_k(i))) x c_tensor x Gradient_s(P_k(j))], where k=x,y,z
   !  Reference: see Eqs. (7) and (8) in:
   !  A.-C. Hladky-Hennion
   !  "Finite element analysis of the propagation of acoustic waves in waveguides,"
   !  Journal of Sound and Vibration, vol. 194, no. 2, pp. 119-136, 1996.
   !
   ! K_ij_sig_tau = c_sig_b_c_tau G_ij_bc
   ! See chapter 9 for details

   do vnd_i=1,P2_NODES_PER_EL         !  v_i nodes
      do v_xyz_i=1,3              !  v_i components  (sig)
         i_p = 3*(vnd_i-1) + v_xyz_i  !  row index

         do und_j=1,P2_NODES_PER_EL         ! u_j nodes
            do u_xyz_j=1,3              ! u_j components  (tau)
               j_p = 3*(und_j-1) + u_xyz_j  ! column index

               ! Do contraction  over stiffness tensor and Gmat
               ! b and c are the inner indices in the contraction as in chap 9

               t_K_ij = C_ZERO
               do b_xyz=1,3
                  vgt_I = S_index(b_xyz,v_xyz_i)  ! S_index is symmetric so order doesn't matter

                  do c_xyz=1,3
                     vgt_J = S_index(c_xyz,u_xyz_j)

                     t_stiffC_IJ = stiffness_C_IJ(vgt_I,vgt_J)

                     t_G_ijbc = Gmat_ij_bc(vnd_i,und_j, b_xyz,c_xyz)

                     ! Account for fact that the u_j_z dof represents the physical (u_{j,z} / i)
                     ! So we need to multiply every u_j_z by i, and every v_j_z^* by -i
                     ! See JLT paper
                     if (v_xyz_i == 3) then
                        t_G_ijbc = -C_IM_ONE* t_G_ijbc
                     endif

                     if (u_xyz_j == 3) then
                        t_G_ijbc = C_IM_ONE* t_G_ijbc
                     endif

                     t_K_ij = t_K_ij + t_G_ijbc * t_stiffC_IJ

                  enddo
               enddo
                mat_K(i_p,j_p) = t_K_ij
            enddo
         enddo
      enddo
   enddo


end

