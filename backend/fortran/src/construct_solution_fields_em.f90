
 !  On entry, evecs_raw is the raw eigenvectors from the arpack soln

 !  On exit:
 !  evecs_raw(*,i) : contains the imaginary and real parts of the solution for points such that cscmat%m_eqs(i) /= 0
 !  sol(i) : contains solution for all points indexed as sol(xyz_comp, 23 nodes per el, n_modes, n_msh_elts)

 !  This is 2D 3-vector component FEM:



 !  Eigenmodes stored in v_evals_beta and xy_ref are reordered according to v_eig_index to sort by largest eigenvalue

#include "numbat_decl.h"

subroutine construct_solution_fields_em (bdy_cdn, shift_ksqr, n_modes, mesh_raw, entities, cscmat, pbcs, &
   n_core, bloch_vec, v_evals_beta, evecs_raw, sol, mode_poln_fracs, nberr)

   use numbatmod
   use class_MeshRaw
   use class_SparseCSC
   use class_PeriodicBCs
   use class_BasisFunctions

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities
   type(SparseCSC) :: cscmat
   type(PeriodicBCs) :: pbcs
   type(BasisFunctions) :: basfuncs

   integer(8) bdy_cdn, n_modes
   integer(8) n_core(2)
   complex(8) shift_ksqr
   double precision bloch_vec(2)

   complex(8) evecs_raw(cscmat%n_dof,n_modes)

   type(NBError), intent(out) :: nberr


   !  sol(3, 1..P2_NODES_PER_EL,n_modes, mesh_raw%n_msh_el)          contains the values of the 3 components at P2 interpolation nodes
   !  sol(3, P2_NODES_PER_EL+1..P2_NODES_PER_EL+7,n_modes, mesh_raw%n_msh_el) contains the values of Ez component at P3 interpolation nodes (per element: 6 edge-nodes and 1 interior node)
   complex(8) sol(3,P2_NODES_PER_EL+7,n_modes,mesh_raw%n_msh_el)
   complex(8) v_evals_beta(n_modes)
   complex(8) mode_poln_fracs(4,n_modes)

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   !  Local variables
   !  32-but integers for BLAS and LAPACK

   integer(8) v_eig_index(n_modes)

   double precision mode_comp(4)
   integer(8) el_nds_i(P2_NODES_PER_EL), phi_vec_map(4,3,N_DDL_T)
   double precision xy_nds_P2(2,P2_NODES_PER_EL), el_nds_xy(2,P2_NODES_PER_EL)
   complex(8) sol_el(3,P2_NODES_PER_EL+7)

   double precision phi_P1_ref(P1_NODES_PER_EL)
   double precision gradt_P1_ref(2,P1_NODES_PER_EL), gradt_P1_act(2,P1_NODES_PER_EL)

   double precision phi_P2_ref(P2_NODES_PER_EL)
   double precision gradt_P2_ref(2,P2_NODES_PER_EL), gradt_P2_act(2,P2_NODES_PER_EL)

   double precision phi_P3_ref(P3_NODES_PER_EL)
   double precision gradt_P3_ref(2,P3_NODES_PER_EL), gradt_P3_act(2,P3_NODES_PER_EL)

   double precision vec_phi_j(2), curl_phi_j, phi_z_j

   double precision xy_ref(2), xy_act(2)
   double precision mat_B(2,2)
   double precision mat_T(2,2)

   complex(8) val_exp(N_ENTITY_PER_EL)

   logical is_curved
   integer(8) j, k, i1, m, nd_i, typ_e
   integer(8) debug, i_sol_max
   integer(8) i_el, md_i, md_i2, jtest, jp, ind_jp, j_eq
   double precision det
   complex(8) z_tmp1, z_tmp2, z_sol_max
   integer(8) nd_lab


   call nberr%reset()
   errco = 0
   debug = 0


   ! Adjust evals to unshifted values and determine reordering
   call rescale_and_sort_eigensolutions(n_modes, shift_ksqr, v_evals_beta, v_eig_index)

   !  Coordinates of the P2 Lagrange interpolation nodes for the unit triangle
   call get_P2_node_locations(P2_NODES_PER_EL, xy_nds_P2)

   mode_poln_fracs  =  D_ZERO

   do md_i=1,n_modes
      md_i2 = v_eig_index(md_i)   !  index of the mode in eigenvalue sorted sequence


      z_sol_max =  D_ZERO         !  value and loc of max field modulus
      i_sol_max = 0

      do i_el=1,mesh_raw%n_msh_el
         typ_e = mesh_raw%el_material(i_el)


         mode_comp =  D_ZERO

         do nd_i=1,P2_NODES_PER_EL
            nd_lab = mesh_raw%elnd_to_mesh(nd_i,i_el)           ! Global label of this node on this elt
            el_nds_i(nd_i) = nd_lab
            el_nds_xy(:,nd_i) = mesh_raw%v_nd_xy(:,nd_lab)      ! xy coords of this node
         enddo

         val_exp =  D_ONE

         if (bdy_cdn == BCS_PERIODIC) then
            call make_pbc_phase_shifts(mesh_raw, entities, pbcs, i_el, bloch_vec, val_exp)
         endif


         call basfuncs%make_phi_vector_map(el_nds_i)

         call make_phi_vector_map (el_nds_i, phi_vec_map)  !  get P2 basis function

         !call is_curved_elem_tri (P2_NODES_PER_EL, el_nds_xy, is_curved, r_tmp1)  !  determine if current element has curved face.
         !Can this ever happen?

         is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, el_nds_xy)

         sol_el = D_ZERO

         do nd_i=1,P2_NODES_PER_EL

            xy_ref = xy_nds_P2(:, nd_i)

            !call basfuncs%evaluate_at_position(i_el, xy_ref, is_curved, el_nds_xy, nberr)
            !RET_ON_NBERR_UNFOLD(nberr)

            !  Elements and gradients for the P1, P2, P3 basis functions
            call phi1_2d_mat (xy_ref, phi_P1_ref, gradt_P1_ref)
            call phi2_2d_mat (xy_ref, phi_P2_ref, gradt_P2_ref)
            call phi3_2d_mat (xy_ref, phi_P3_ref, gradt_P3_ref)

            if (.not. is_curved ) then
               !  Rectilinear element
               call jacobian_p1_2d (xy_ref, el_nds_xy, P2_NODES_PER_EL, xy_act, det, mat_B, mat_T, errco, emsg)
               RETONERROR(errco)

               if (det <= 0 .and. debug == 1) then
                  write(*,*) "   !!!"
                  write(*,*) "array_sol: det <= 0: i_el, det ", i_el, det
               endif

            else
               !  Isoparametric element, 2024-06-14 fix
               call jacobian_p2_2d (el_nds_xy, P2_NODES_PER_EL, phi_P2_ref, gradt_P2_ref, &
                  xy_act, det, mat_B, mat_T, errco, emsg)
               RETONERROR(errco)
            endif


            !  grad_i  = gradient on the actual triangle
            !  grad_i  = Transpose(mat_T)*grad_i0
            !  Calculation of the matrix-matrix product:
            call DGEMM('Transpose','N', 2, 3,  2, D_ONE, mat_T, 2, &
               gradt_P1_ref, 2, D_ZERO, gradt_P1_act, 2)
            call DGEMM('Transpose','N', 2, 6,  2, D_ONE, mat_T, 2, &
               gradt_P2_ref, 2, D_ZERO, gradt_P2_act, 2)
            call DGEMM('Transpose','N', 2, 10, 2, D_ONE, mat_T, 2, &
               gradt_P3_ref, 2, D_ZERO, gradt_P3_act, 2)

            !  Contribution to the transverse component
            do jtest=1,N_DDL_T
               do j_eq=1,3
                  jp = entities%v_tags(jtest,i_el)
                  ind_jp = cscmat%m_eqs(j_eq,jp)
                  if (ind_jp > 0) then
                     m  = phi_vec_map(2, j_eq, jtest)
                     if (m == nd_i) then

                        !  nd_i correspond to a P2 interpolation node
                        !  The contribution is nonzero only when m=nd_i.
                        !  Determine the basis vector

                        call make_phi_vector_basis(j_eq, jtest, phi_vec_map, phi_P2_ref, &
                           gradt_P1_act, gradt_P2_act, vec_phi_j, curl_phi_j)
                        z_tmp1 = evecs_raw(ind_jp, md_i2)* val_exp(jtest)


                        do j=1,2
                           z_tmp2 = z_tmp1 * vec_phi_j(j)
                           sol_el(j,nd_i) = sol_el(j,nd_i) + z_tmp2

                           if (m /= nd_i .and. abs(z_tmp2) > 1.0d-7) then
                              write(*,*)
                              write(*,*) i_el, nd_i, m, abs(z_tmp2)
                              write(*,*) "vec_phi_j = ", vec_phi_j
                              write(*,*) "xy_ref = ", xy_ref
                              write(*,*) "xy_nds_P2 = ", (xy_nds_P2(k,nd_i),k=1,2)
                              write(*,*) "phi_P2_ref = ", phi_P2_ref
                           endif

                        enddo

                     endif
                  endif
               enddo
            enddo

            !  Contribution to the longitudinal component
            !  The initial P3 value of Ez isinterpolated over P2 nodes
            do jtest=N_DDL_T+1,N_ENTITY_PER_EL

               do j_eq=1,1
                  jp = entities%v_tags(jtest,i_el)
                  ind_jp = cscmat%m_eqs(j_eq,jp)
                  if (ind_jp > 0) then

                     m  = jtest-N_DDL_T
                     phi_z_j = phi_P3_ref(m)


                     z_tmp2 = evecs_raw(ind_jp, md_i2) * val_exp(jtest) * phi_z_j
                     sol_el(3,nd_i) = sol_el(3,nd_i) + z_tmp2
                  endif
               enddo
            enddo

            do j=1,3
               z_tmp2 = sol_el(j,nd_i)
               sol(j,nd_i,md_i,i_el) = z_tmp2
               if (abs(z_sol_max) < abs(z_tmp2)) then  !  found a new max (by component not total?)
                  z_sol_max = z_tmp2
                  i_sol_max = mesh_raw%elnd_to_mesh(nd_i,i_el)
               endif
            enddo

            !  Contribution of the element i_el to the mode component
            mode_comp(1:3) = mode_comp(1:3) + abs(sol_el(1:3,nd_i))**2

         enddo

         !  Saving the P3 values of Ez at: the 6 edge nodes and the interior node
         do nd_i=P2_NODES_PER_EL+1,P2_NODES_PER_EL+7

            sol_el(1:3,nd_i) = D_ZERO

            jtest = N_DDL_T+nd_i-P2_NODES_PER_EL+3
            j_eq = 1
            jp = entities%v_tags(jtest,i_el)
            ind_jp = cscmat%m_eqs(j_eq,jp)

            if (ind_jp > 0) then
               sol_el(3,nd_i) = evecs_raw(ind_jp, md_i2)* val_exp(jtest)
            endif

            sol(1:3,nd_i,md_i,i_el) =  sol_el(1:3,nd_i)

         enddo

         !  Avarage values

         mode_comp(1:3) = mode_comp(1:3) * abs(det)/dble(P2_NODES_PER_EL)

         !  Add the contribution of the element i_el to the mode component
         mode_poln_fracs(1:3, md_i) = mode_poln_fracs(1:3, md_i) + mode_comp(1:3)

         if (typ_e == n_core(1) .or. typ_e == n_core(2)) then
            mode_poln_fracs(4,md_i) = mode_poln_fracs(4,md_i) + mode_comp(1) + mode_comp(2) + mode_comp(3)
         endif

      enddo

      !  Total energy and normalization
      z_tmp2 = mode_poln_fracs(1,md_i) + mode_poln_fracs(2,md_i) + mode_poln_fracs(3,md_i)
      !if (abs(z_tmp2) < 1.0d-10) then
      if (abs(z_tmp2) < 1.0d-20) then ! 11/12/2024, trying to allow thin triangle element
         write(*,*) "array_sol: the total energy ",        "is too small : ", z_tmp2
         write(*,*) "array_sol: md_i md_i2 = ", md_i, md_i2
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif

      mode_poln_fracs(:,md_i) = mode_poln_fracs(:,md_i) / z_tmp2

      !  Check if the eigenvector is nonzero
      !if (abs(z_sol_max) < 1.0d-10) then
      if (abs(z_sol_max) < 1.0d-20) then ! 11/12/2024, trying to allow thin triangle element
         z_sol_max = z_tmp2
         write(*,*) "array_sol: z_sol_max is too small"
         write(*,*) "array_sol: z_sol_max = ", z_sol_max
         write(*,*) "md_i, md_i2, n_modes = ", md_i, md_i2, n_modes
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif

      !  Normalization so that the maximum fi_eld component is 1
      do i_el=1,mesh_raw%n_msh_el
         do nd_i=1,P2_NODES_PER_EL
            i1 = mesh_raw%elnd_to_mesh(nd_i,i_el)

            sol(1:3,nd_i,md_i,i_el) = sol(1:3,nd_i,md_i,i_el)/z_sol_max

            i1 = mesh_raw%elnd_to_mesh(nd_i,i_el)
            if (i1 == i_sol_max .and. debug == 1) then
               write(*,*) "array_sol:"
               write(*,*) "md_i, i1, i_el = ", md_i, i1, i_el
               write(*,*) "array_sol: Fi_eld normalisaion point:"
               write(*,*) "x = ", dble(mesh_raw%v_nd_xy(1,i1))
               write(*,*) "y = ", dble(mesh_raw%v_nd_xy(2,i1))
               write(*,*) "i_sol_max = ", i_sol_max
               write(*,*) md_i, i1, i_el, (dble(sol(j,nd_i,md_i,i_el)),j=1,3)
               write(*,*) md_i, i1, i_el, (imag(sol(j,nd_i,md_i,i_el)),j=1,3)
            endif
         enddo


         do nd_i=P2_NODES_PER_EL+1,P2_NODES_PER_EL+7
            sol(1:3,nd_i,md_i,i_el) = sol(1:3,nd_i,md_i,i_el)/z_sol_max
         enddo

      enddo



      evecs_raw(1:cscmat%n_dof,md_i2) = evecs_raw(1:cscmat%n_dof,md_i2)/z_sol_max
   enddo


   !  The z-component must be multiplied by -i*beta to recover the physical,
   !  un-normalised z-component (see Eq. (25) of the JOSAA 2012 paper)
   do md_i=1,n_modes
      sol(3,:,md_i,:) = C_IM_ONE * v_evals_beta(md_i) * sol(3,:,md_i,:)
   enddo


end




subroutine zvec_reorder_by_index(v_src_dest, v_eig_index, num_elts)


   integer(8) :: num_elts
   complex(8) :: v_src_dest(num_elts)
   integer(8) :: v_eig_index(num_elts)

   complex(8) ::  v_tmp(num_elts)
   integer(8) :: j, j1

   do j=1,num_elts
      j1=v_eig_index(j)
      v_tmp(j) = v_src_dest(j1)
   enddo

   v_src_dest = v_tmp

end subroutine

subroutine rescale_and_sort_eigensolutions(n_modes, shift_ksqr, v_evals_beta, v_eig_index)

   integer(8), intent(in) :: n_modes
   complex(8), intent(in) :: shift_ksqr
   complex(8) :: v_evals_beta(n_modes)
   integer(8) :: v_eig_index(n_modes)

   integer(8) i

   complex(8) z_beta

   !TODO: make a function. Turn beta^2 raw eig into actual beta
   do i=1,n_modes
      !  z_tmp0 = v_evals_beta(i)
      !  z_tmp = 1.0d0/z_tmp0+shift_ksqr
      !  z_beta = sqrt(z_tmp)

      z_beta = sqrt(1.0d0/v_evals_beta(i)+shift_ksqr )

      !  Mode classification - we want the forward propagating mode
      if (abs(imag(z_beta)/z_beta) .lt. 1.0d-8) then
         !  re(z_beta) > 0 for forward propagating mode
         if (dble(z_beta) .lt. 0) z_beta = -z_beta
      else
         !  im(z_beta) > 0 for forward decaying evanescent mode  !rarely relevant for us
         if (imag(z_beta) .lt. 0) z_beta = -z_beta
      endif

      v_evals_beta(i) = z_beta
   enddo

   !  order v_evals_beta by magnitudes and store in v_eig_index
   call z_indexx (n_modes, v_evals_beta, v_eig_index)

   ! Apply the reordering
   call zvec_reorder_by_index(v_evals_beta, v_eig_index, n_modes)

end subroutine


subroutine make_pbc_phase_shifts(mesh_raw, entities, pbcs, i_el, bloch_vec, val_exp)

   use numbatmod
   use class_MeshRaw
   use class_PeriodicBCs

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities
   type(PeriodicBCs) :: pbcs
   double precision bloch_vec(2)

   complex(8) val_exp(N_ENTITY_PER_EL)

   integer(8) i_el, nd_i, nd_lab, j1, jp, j, k

   integer(8) el_nds_i(P2_NODES_PER_EL)



   complex(8) r_tmp1
   double precision delta_xy_ref(2)
   double precision ddot

   val_exp = D_ONE

   do nd_i=1,P2_NODES_PER_EL
      nd_lab = mesh_raw%elnd_to_mesh(nd_i,i_el)
      k = pbcs%iperiod_N(nd_lab)
      if (k /= 0) nd_lab=k
      el_nds_i(nd_i) = nd_lab
   enddo

   !  val_exp: Bloch mod ephase factor between the origin point and destination point
   !  For a pair of periodic points, one is chosen as origin and the other is the destination
   do j=1,N_ENTITY_PER_EL
      jp = entities%v_tags(j,i_el)
      j1 = pbcs%iperiod_N_E_F(jp)
      if (j1 /= 0) then
         !do k=1,dim_32
         !  delta_xy_ref(k) = entities.v_nd_xy(k,jp) - entities.v_nd_xy(k,j1)
         !enddo
         delta_xy_ref = entities%v_xy(:,jp) - entities%v_xy(:,j1)
         r_tmp1 = ddot(2, bloch_vec, 1, delta_xy_ref, 1)
         val_exp(j) = exp(C_IM_ONE * r_tmp1)
      endif
   enddo

end subroutine

