
!   On exit:
!     sol_0(*,i) : contains the imaginary and real parts of the solution for points such that ineq(i) /= 0
!     sol(i) : contains solution for all points

!     This is 2D 3-vector component FEM:
!        The dimension of the geometric domain is : dim_32 = 2
!        The dimension of the vector field is : dim2 = 3
!

!    Eigenmodes stored in v_eigs_beta and XX are reordered according to iindex to sort by largest eigenvalue

subroutine array_sol (i_cond, num_modes, n_msh_el, n_msh_pts, &
   n_ddl, neq, nnodes, n_core, bloch_vec, iindex, table_nod, &
   table_N_E_F, type_el, ineq,      ip_period_N, ip_period_N_E_F, &
   mesh_xy, x_N_E_F, v_eigs_beta, mode_pol, sol_0, sol)

   use numbatmod

   integer*8 i_cond, num_modes, n_msh_el, n_msh_pts, n_ddl
   integer*8 neq, nnodes
   integer*8 n_core(2)
   integer*8 type_el(n_msh_el)
   integer*8 ineq(3,n_ddl)   ! bc info
   integer*8 iindex(*)
   integer*8 ip_period_N(n_msh_pts), ip_period_N_E_F(n_ddl)
   integer*8 table_nod(nnodes,n_msh_el), table_N_E_F(14,n_msh_el)
   double precision bloch_vec(2), mesh_xy(2,n_msh_pts)
   double precision x_N_E_F(2,n_ddl)
   complex*16 sol_0(neq,num_modes)

!     sol(3, 1..nnodes,num_modes, n_msh_el)          contains the values of the 3 components at P2 interpolation nodes
!     sol(3, nnodes+1..nnodes+7,num_modes, n_msh_el) contains the values of Ez component at P3 interpolation nodes (per element: 6 edge-nodes and 1 interior node)
   complex*16 sol(3,nnodes+7,num_modes,n_msh_el)
   complex*16 v_eigs_beta(num_modes)
   complex*16 mode_pol(4,num_modes)


!     Local variables
!      integer*8 nnodes_0, nddl_0, nddl_t
!     32-but integers for BLAS and LAPACK
   integer*8 nddl_0
   integer*4 dim_32
!      parameter (nnodes_0 = 6)
   parameter (nddl_0 = 14)
!      parameter (nddl_t=4)
   parameter (dim_32=2)
!
   double precision mode_comp(4)
   integer*8 nod_el_p(nnodes_0), basis_list(4,3,nddl_t)
   double precision xn(dim_32,nnodes_0), el_xy(dim_32,nnodes_0)
   complex*16 sol_el(3,nnodes_0+7)

   double precision phi1_list(3), grad1_mat0(dim_32,3)
   double precision grad1_mat(dim_32,3)

   double precision phi2_list(6), grad2_mat0(dim_32,6)
   double precision grad2_mat(dim_32,6)

   double precision phi3_list(10), grad3_mat0(dim_32,10)
   double precision grad3_mat(dim_32,10)

   double precision vec_phi_j(dim_32), curl_phi_j, phi_z_j

   double precision xx(dim_32), xx_g(dim_32)
   double precision delta_xx(dim_32)
   double precision mat_B(dim_32,dim_32)
   double precision mat_T(dim_32,dim_32)

   double complex val_exp(nddl_0)

   integer*8 is_curved
   integer*8 j, k, i1, j1, m, inod, typ_e
   integer*8 debug, i_sol_max
   integer*8 iel, i_mode, i_mode2, jtest, jp, ind_jp, j_eq
   double precision ddot, det, r_tmp1
   complex*16 z_tmp1, z_tmp2, z_sol_max
   debug = 0


!   Reorder eigenvectors by iindex
   call zvec_reorder_by_index(v_eigs_beta, iindex, num_modes)


   !do j=1,num_modes
   !   j1=iindex(j)
   !   v_tmp(j) = v_eigs_beta(j1)
   !enddo
   !v_eigs_beta = v_tmp

!     Coordinates of the P2 Lagrange interpolation nodes for the unit triangle
   call interp_nod_2d (nnodes, xn)

   do i_mode=1,num_modes
      i_mode2 = iindex(i_mode)   ! index of the next largest eigenvalue

      mode_pol(:,i_mode) =  D_ZERO

      z_sol_max =  D_ZERO   ! value and loc of max field modulus
      i_sol_max = 0

      do iel=1,n_msh_el
         typ_e = type_el(iel)


         mode_comp =  D_ZERO

         do inod=1,nnodes
            j = table_nod(inod,iel)
            nod_el_p(inod) = j
            el_xy(:,inod) = mesh_xy(:,j)
         enddo

         val_exp =  D_ONE

         if (i_cond == BCS_PERIODIC) then
!           Periodic boundary condition
            do inod=1,nnodes
               j = table_nod(inod,iel)
               k = ip_period_N(j)
               if (k /= 0) j=k
               nod_el_p(inod) = j
            enddo

!           val_exp: Bloch mod ephase factor between the origin point and destination point
!           For a pair of periodic points, one is chosen as origin and the other is the destination
            do j=1,nddl_0
               jp = table_N_E_F(j,iel)
               j1 = ip_period_N_E_F(jp)
               if (j1 /= 0) then
                  !do k=1,dim_32
                  !   delta_xx(k) = x_N_E_F(k,jp) - x_N_E_F(k,j1)
                  !enddo
                  delta_xx(:) = x_N_E_F(:,jp) - x_N_E_F(:,j1)
                  r_tmp1 = ddot(dim_32, bloch_vec, 1, delta_xx, 1)
                  val_exp(j) = exp(C_IM_ONE * r_tmp1)
               endif
            enddo
         endif  ! BCS_PERIODIC

         call basis_ls (nod_el_p, basis_list)  ! get P2 basis function

         call curved_elem_tri (nnodes, el_xy, is_curved, r_tmp1)  ! determine if current element has curved face. Can this ever happen?
!         P2 Lagrange Interpolation nodes for the unit triangle
!         xn   = coordinate on the reference triangle
!          do inod=1,nnodes+7
!            do j=1,3
!              sol_el(j,inod) = 0.00
!            enddo
!          enddo

         do inod=1,nnodes

            xx =xn(:, inod)
            sol_el(:, inod) = D_ZERO

            ! Elements and gradients for the P1, P2, P3 basis functions
            call phi1_2d_mat (xx, phi1_list, grad1_mat0)
            call phi2_2d_mat (xx, phi2_list, grad2_mat0)
            call phi3_2d_mat (xx, phi3_list, grad3_mat0)

            if (is_curved == 0) then
!             Rectilinear element
               call jacobian_p1_2d (xx, el_xy, nnodes, xx_g, det, mat_B, mat_T)

               if (det <= 0 .and. debug == 1) then
                  write(*,*) "   !!!"
                  write(*,*) "array_sol: det <= 0: iel, det ", iel, det
               endif

            else
!             Isoparametric element, 2024-06-14 fix
               call jacobian_p2_2d (el_xy, nnodes, phi2_list, grad2_mat0, xx_g, det, mat_B, mat_T)
            endif

!             if(abs(det) < 1.0d-10) then
            if(abs(det) < 1.0d-20) then
               write(*,*)
               write(*,*) "   ???"
               write(*,*) "array_sol: det = 0 : iel, det = ", iel, det
               write(*,*) "array_sol: Aborting..."
               stop
            endif

!           grad_i  = gradient on the actual triangle
!           grad_i  = Transpose(mat_T)*grad_i0
!           Calculation of the matrix-matrix product:
            call DGEMM('Transpose','N', dim_32, 3,  dim_32, D_ONE, mat_T, dim_32, grad1_mat0, dim_32, D_ZERO, grad1_mat, dim_32)
            call DGEMM('Transpose','N', dim_32, 6,  dim_32, D_ONE, mat_T, dim_32, grad2_mat0, dim_32, D_ZERO, grad2_mat, dim_32)
            call DGEMM('Transpose','N', dim_32, 10, dim_32, D_ONE, mat_T, dim_32, grad3_mat0, dim_32, D_ZERO, grad3_mat, dim_32)
!
!           Contribution to the transverse component
            do jtest=1,nddl_t
               do j_eq=1,3
                  jp = table_N_E_F(jtest,iel)
                  ind_jp = ineq(j_eq,jp)
                  if (ind_jp > 0) then
                     m  = basis_list(2, j_eq, jtest)
                     if (m == inod) then
!                   !  inod correspond to a P2 interpolation node
!                                           The contribution is nonzero only when m=inod.
!                 Determine the basis vector

                        call basis_vec (j_eq, jtest, basis_list, phi2_list, grad1_mat, grad2_mat, vec_phi_j, curl_phi_j)
                        z_tmp1 = sol_0(ind_jp, i_mode2)* val_exp(jtest)


                        do j=1,dim_32
                           z_tmp2 = z_tmp1 * vec_phi_j(j)
                           sol_el(j,inod) = sol_el(j,inod) + z_tmp2

                           if (m /= inod .and. abs(z_tmp2) > 1.0d-7) then
                              write(*,*)
                              write(*,*) iel, inod, m, abs(z_tmp2)
                              write(*,*) "vec_phi_j = ", vec_phi_j
                              write(*,*) "xx = ", xx
                              write(*,*) "xn = ", (xn(k,inod),k=1,dim_32)
                              write(*,*) "phi2_list = ", phi2_list
                           endif

                        enddo

                     endif
                  endif
               enddo
            enddo

!           Contribution to the longitudinal component
!             !  The initial P3 value of Ez isinterpolated over P2 nodes
            do jtest=nddl_t+1,nddl_0
!               ! 3
               do j_eq=1,1
                  jp = table_N_E_F(jtest,iel)
                  ind_jp = ineq(j_eq,jp)
                  if (ind_jp > 0) then
                     !z_tmp1 = sol_0(ind_jp, i_mode2)
                     m  = jtest-nddl_t
                     phi_z_j = phi3_list(m)
                     !!z_tmp1 = z_tmp1 * val_exp(jtest)
                     !z_tmp2 = z_tmp1 * phi_z_j

                     z_tmp2 = sol_0(ind_jp, i_mode2) * val_exp(jtest) * phi_z_j
                     sol_el(3,inod) = sol_el(3,inod) + z_tmp2
                  endif
               enddo
            enddo

            do j=1,3
               z_tmp2 = sol_el(j,inod)
               sol(j,inod,i_mode,iel) = z_tmp2
               if (abs(z_sol_max) < abs(z_tmp2)) then  ! found a new max (by component not total?)
                  z_sol_max = z_tmp2
                  i_sol_max = table_nod(inod,iel)
               endif
            enddo

!           Contribution of the element iel to the mode component
            !do j=1,3
            !!   mode_comp(j) = mode_comp(j) + abs(sol_el(j,inod))**2
            !enddo
            mode_comp(1:3) = mode_comp(1:3) + abs(sol_el(1:3,inod))**2

         enddo
!ccccccccc
!         Saving the P3 values of Ez at: the 6 edge nodes and the interior node
         do inod=nnodes+1,nnodes+7
            !do j=1,3
            !   sol_el(j,inod) = 0.00
            !enddo
            sol_el(1:3,inod) = D_ZERO

            jtest = nddl_t+inod-nnodes+3
            j_eq = 1
            jp = table_N_E_F(jtest,iel)
            ind_jp = ineq(j_eq,jp)

            if (ind_jp > 0) then
               !z_tmp1 = sol_0(ind_jp, i_mode2)
               !z_tmp1 = z_tmp1 * val_exp(jtest)
               !sol_el(3,inod) = z_tmp1
               sol_el(3,inod) = sol_0(ind_jp, i_mode2)* val_exp(jtest)
            endif

            !do j=1,3
            !!   z_tmp2 = sol_el(j,inod)
             !  sol(j,inod,i_mode,iel) = z_tmp2
            !enddo
            sol(1:3,inod,i_mode,iel) =  sol_el(1:3,inod)

         enddo

!         Avarage values
         !do j=1,3
         !   mode_comp(j) = abs(det)*mode_comp(j)/dble(nnodes)
         !enddo
         mode_comp(1:3) = mode_comp(1:3) * abs(det)/dble(nnodes)

!         Add the contribution of the element iel to the mode component
         !do j=1,3
         !   mode_pol(j,i_mode) = mode_pol(j,i_mode) + mode_comp(j)
         !enddo
         mode_pol(1:3, i_mode) = mode_pol(1:3, i_mode) + mode_comp(1:3)

         if (typ_e == n_core(1) .or. typ_e == n_core(2)) then
            !z_tmp2 = mode_comp(1) + mode_comp(2)         + mode_comp(3)
            !mode_pol(4,i_mode) = mode_pol(4,i_mode) + z_tmp2
            mode_pol(4,i_mode) = mode_pol(4,i_mode) + mode_comp(1) + mode_comp(2) + mode_comp(3)
         endif

      enddo

!       Total energy and normalization
      z_tmp2 = mode_pol(1,i_mode) + mode_pol(2,i_mode) + mode_pol(3,i_mode)
      if (abs(z_tmp2) < 1.0d-10) then
         write(*,*) "array_sol: the total energy ",        "is too small : ", z_tmp2
         write(*,*) "array_sol: i_mode i_mode2 = ", i_mode, i_mode2
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif

      ! do j=1,3
      !    mode_pol(j,i_mode) = mode_pol(j,i_mode) / z_tmp2
      ! enddo

      ! j=4
      ! mode_pol(j,i_mode) = mode_pol(j,i_mode) / z_tmp2
      mode_pol(:,i_mode) = mode_pol(:,i_mode) / z_tmp2

!       Check if the eigenvector is nonzero
      if (abs(z_sol_max) < 1.0d-10) then
         z_sol_max = z_tmp2
         write(*,*) "array_sol: z_sol_max is too small"
         write(*,*) "array_sol: z_sol_max = ", z_sol_max
         write(*,*) "i_mode, i_mode2, num_modes = ", i_mode, i_mode2, num_modes
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif

!       Normalization so that the maximum field component is 1
      do iel=1,n_msh_el
         do inod=1,nnodes
            i1 = table_nod(inod,iel)

            !do j=1,3
            !!   z_tmp1 = sol(j,inod,i_mode,iel)/z_sol_max
             !  sol(j,inod,i_mode,iel) = z_tmp1
            !enddo

            sol(1:3,inod,i_mode,iel) = sol(1:3,inod,i_mode,iel)/z_sol_max

            i1 = table_nod(inod,iel)
            if (i1 == i_sol_max .and. debug == 1) then
               write(*,*) "array_sol:"
               write(*,*) "i_mode, i1, iel = ", i_mode, i1, iel
               write(*,*) "array_sol: Field normalisaion point:"
               write(*,*) "x = ", dble(mesh_xy(1,i1))
               write(*,*) "y = ", dble(mesh_xy(2,i1))
               write(*,*) "i_sol_max = ", i_sol_max
               write(*,*) i_mode, i1, iel, (dble(sol(j,inod,i_mode,iel)),j=1,3)
               write(*,*) i_mode, i1, iel, (imag(sol(j,inod,i_mode,iel)),j=1,3)
            endif
         enddo

!ccccccccc
         do inod=nnodes+1,nnodes+7
            !do j=1,3
            !   z_tmp1 = sol(j,inod,i_mode,iel)/z_sol_max
            !   sol(j,inod,i_mode,iel) = z_tmp1
            !enddo
            sol(1:3,inod,i_mode,iel) = sol(1:3,inod,i_mode,iel)/z_sol_max
         enddo

!ccccccccc
      enddo

      !do j=1,neq
      !   z_tmp1 = sol_0(j,i_mode2)/z_sol_max
      !   sol_0(j,i_mode2) = z_tmp1
      !enddo

      sol_0(1:neq,i_mode2) = sol_0(1:neq,i_mode2)/z_sol_max
   enddo
!
   return
end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zvec_reorder_by_index(v_src_dest, iindex, num_elts)


   integer*8 :: num_elts
   complex*16 :: v_src_dest(num_elts)
   integer*8 :: iindex(num_elts)

   complex*16 ::  v_tmp(num_elts)
   integer*8 :: j, j1

   do j=1,num_elts
      j1=iindex(j)
      v_tmp(j) = v_src_dest(j1)
   enddo

   v_src_dest = v_tmp

end subroutine