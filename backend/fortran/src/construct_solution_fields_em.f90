
 !  On entry, evecs_raw is the raw eigenvectors from the arpack eigensolution

 !  On exit:
 !  evecs_raw(*,i) : contains the imaginary and real parts of the arpack eigensolution for all active dof
 !  evecs_final(i) : contains eigensolution for all field points indexed as evecs_final(xyz_comp, 23 nodes per el, n_modes, n_msh_eltsts)

 !  Eigenmodes stored in v_evals_beta and xy_ref are reordered according to v_eig_index to sort by largest eigenvalue

#include "numbat_decl.h"

subroutine construct_solution_fields_em (bdy_cdn, shift_ksqr, n_modes, mesh, entities, &
   cscmat, pbcs, bloch_vec, v_evals_beta, evecs_raw, evecs_final, mode_poln_fracs, nberr)

   use numbatmod
   use class_Mesh
   use class_SparseCSC_EM
   use class_PeriodicBCs
   use class_BasisFunctions

   type(MeshEM) :: mesh
   type(MeshEntities) :: entities
   type(SparseCSC_EM) :: cscmat
   type(PeriodicBCs) :: pbcs

   integer(8) bdy_cdn, n_modes
   complex(8) shift_ksqr
   double precision bloch_vec(2)

   complex(8) evecs_raw(cscmat%n_dof,n_modes)

   type(NBError), intent(out) :: nberr




   complex(8) evecs_final(3,N_DOF_PER_EL,n_modes,mesh%n_msh_elts)
   complex(8) v_evals_beta(n_modes)
   complex(8) mode_poln_fracs(4,n_modes)

   integer(8) md_i, md_i_sorted


   !  Local variables

   integer(8) v_eig_index(n_modes)

   ! Adjust evals to unshifted values and determine reordering
   call rescale_and_sort_eigensolutions(n_modes, shift_ksqr, v_evals_beta, v_eig_index)

   mode_poln_fracs  =  D_ZERO

   do md_i=1,n_modes
      md_i_sorted = v_eig_index(md_i)   !  index of the mode in eigenvalue sorted sequence

      call construct_em_mode_i(md_i, md_i_sorted, n_modes,  mesh, entities, cscmat, &
         evecs_raw, evecs_final, mode_poln_fracs, nberr)

      RET_ON_NBERR(nberr)

   enddo


   !  The z-component must be multiplied by -i*beta to recover the physical,
   !  un-normalised z-component (see Eq. (25) of the JOSAA 2012 paper)
   do md_i=1,n_modes
      evecs_final(3,:,md_i,:) = C_IM_ONE * v_evals_beta(md_i) * evecs_final(3,:,md_i,:)
   enddo

end


 !  evecs_final(3, 1..P2_NODES_PER_EL,n_modes, mesh%n_msh_elts)          contains the values of the 3 components at P2 interpolation nodes
 !  evecs_final(3, P2_NODES_PER_EL+1..N_DOF_PER_EL,n_modes, mesh%n_msh_elts) contains the values of Ez component at P3 interpolation nodes (per element: 6 edge-nodes and 1 interior node)

subroutine construct_em_mode_i(md_i, md_i_sorted, n_modes, mesh, entities, cscmat, &
   evecs_raw, evecs_final, mode_poln_fracs, nberr)
   use numbatmod
   use class_Mesh
   use class_SparseCSC_EM
   use class_BasisFunctions


   integer(8) md_i, md_i_sorted, n_modes

   type(MeshEM) :: mesh
   type(MeshEntities) :: entities
   type(SparseCSC_EM) :: cscmat

   complex(8) evecs_raw(cscmat%n_dof,n_modes)


   complex(8) evecs_final(3,N_DOF_PER_EL,n_modes,mesh%n_msh_elts)
   complex(8) mode_poln_fracs(4,n_modes)

   type(NBError), intent(out) :: nberr

   !----------------

   type(BasisFunctions) :: basfuncs
   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision xy_nds_P2(2,P2_NODES_PER_EL), el_nds_xy(2,P2_NODES_PER_EL)

   double precision mode_comp(4)

   complex(8) evecs_final_el(3,N_DOF_PER_EL) ! evecs_final solution for one mode and elt

   double precision vec_phi_x(2), curlt_phi_x, phi_P3_x   !TODO: move these into basfuncs

   double precision xy_ref(2)
   logical is_curved
   integer(8) xyz_i, nd_p3, nd_p2
   integer(8) i_el,  ety_j, ety_tag, absdof, locdof_j
   integer(8) owning_node, ety_p3
   complex(8) z_tmp2, z_evecs_final_max, E_xyz

   character(len=EMSG_LENGTH) emsg


   !------------------------

   !  Coordinates of the P2 Lagrange interpolation nodes for the unit triangle
   call get_P2_node_locations(xy_nds_P2)

   !
   z_evecs_final_max = D_ZERO   !  value and loc of max field modulus
   !i_evecs_final_max = 0

   do i_el=1,mesh%n_msh_elts

      mode_comp = D_ZERO
      !val_exp = D_ONE
      evecs_final_el = D_ZERO


      call mesh%find_nodes_for_elt(i_el, el_nds_i, el_nds_xy, is_curved)


      ! if (bdy_cdn == BCS_PERIODIC) then
      !    call make_pbc_phase_shifts(mesh, entities, pbcs, i_el, bloch_vec, val_exp)
      ! endif


      call basfuncs%build_vector_elt_map(el_nds_i)


      ! fill evecs_final_el(1:3, 1..P2_NODES)
      do nd_p2=1,P2_NODES_PER_EL
         xy_ref = xy_nds_P2(:, nd_p2)

         !TODO:  create a basfuncs%evaluate_at_P2_node function and remove xy_nds_P2 and get_P2_node_locations
         !       these are used nowhere else

         call basfuncs%evaluate_at_position(i_el, xy_ref, is_curved, el_nds_xy, nberr)
         RET_ON_NBERR(nberr)

         ! transverse part:  evecs_final_el(1:2, 1..P2_NODES)
         do ety_j=1,N_ETY_TRANSVERSE  ! for the transverse field entities on this elt
            ety_tag = entities%v_tags(ety_j, i_el)    ! find the global ety id

            do locdof_j=1,3                            ! the entity can have up to 3 dof
               absdof = cscmat%m_global_dofs(locdof_j, ety_tag)   ! eq num for this dof

               if (absdof > 0) then

                  ! The vector elements are built from P2 scalar functions which are nonzero
                  !  at a P2 node, only for the function corresponding to that node
                  ! So we should only evaluate basis functions which are relevant at that node
                  !TODO: create a basfuncs%get_scalar_index_of_vector_elt function

                  owning_node  = basfuncs%vector_elt_map(2, locdof_j, ety_j)

                  if (owning_node == nd_p2) then

                     call basfuncs%evaluate_vector_elts(locdof_j, ety_j, vec_phi_x, curlt_phi_x)

                     ! pbc version
                     !evecs_final_el(1:2,nd_p2) = evecs_final_el(1:2,nd_p2) + evecs_raw(absdofs, md_i_sorted) * vec_phi_x* val_exp(ety_j)
                     evecs_final_el(1:2,nd_p2) = evecs_final_el(1:2,nd_p2) + evecs_raw(absdof, md_i_sorted) * vec_phi_x

                  endif
               endif
            enddo
         enddo

         !  Longitudinal part:  evecs_final_el(3, 1..P2_NODES)
         !  Contribution to the longitudinal component
         !  The initial P3 value of Ez is interpolated over P2 nodes
         do nd_p3=1, P3_NODES_PER_EL  ! from 5..14 = 10 P3 nodes
            ety_p3 = ETY_TAG_OFFSET_P3_VERTICES + nd_p3

            !do ety_j=N_ETY_TRANSVERSE+1,N_ENTITY_PER_EL  ! from 5..14 = 10 P3 nodes

            ety_tag = entities%v_tags(ety_p3, i_el)

            locdof_j=1                                  ! There is only 1 DOF for each of the P3 nodes
            absdof = cscmat%m_global_dofs(locdof_j, ety_tag)     ! Find the index of this dof at this entity

            if (absdof > 0) then

               !      m  = ety_j-N_ETY_TRANSVERSE
               !     phi_P3_x = basfuncs%phi_P3_ref(m)
               phi_P3_x = basfuncs%phi_P3_ref(nd_p3)

               !pbc version
               !evecs_final_el(3,nd_p2) = evecs_final_el(3,nd_p2) + evecs_raw(absdof, md_i_sorted) * phi_P3_x * val_exp(ety_j)

               evecs_final_el(3,nd_p2) = evecs_final_el(3,nd_p2) + evecs_raw(absdof, md_i_sorted) * phi_P3_x

            endif
         enddo


         ! check if we have a new maximum sized component
         do xyz_i=1,3
            E_xyz = evecs_final_el(xyz_i,nd_p2)
            if (abs(z_evecs_final_max) < abs(E_xyz)) then  !  found a new max
               z_evecs_final_max = E_xyz
               ! i_evecs_final_max = mesh%m_elnd_to_mshpt(nd_p2,i_el)
            endif
         enddo

         !  Contribution of the element i_el to the mode component
         mode_comp(1:3) = mode_comp(1:3) + abs(evecs_final_el(1:3,nd_p2))**2

      enddo  ! end of current P2 node




      !  Average values
      mode_comp(1:3) = mode_comp(1:3) * abs(basfuncs%det)/dble(P2_NODES_PER_EL)


      !  Longitudinal part:  evecs_final_el(3, P3_NODES...)
      !   x and comps of the P3_NODES are left empty
      !  Saving the P3 values of Ez at: the 6 edge nodes and the interior node
      ! This is a straight copy
!      do nd_p3=P2_NODES_PER_EL+1,N_DOF_PER_EL  ! 7..13 = 7 P3 nodes

      do nd_p3=ETY_TAG_OFFSET_P3_NODES,N_DOF_PER_EL  ! 7..13 = 7 P3 nodes


         !evecs_final_el(1:3,nd_p3) = D_ZERO

         !ety_j = N_ETY_TRANSVERSE+nd_p3-P2_NODES_PER_EL+3
         ety_j = nd_p3 + 1  ! make space for the face element
         locdof_j = 1
         ety_tag = entities%v_tags(ety_j,i_el)
         absdof = cscmat%m_global_dofs(locdof_j,ety_tag)

         if (absdof > 0) then
            !evecs_final_el(3,nd_p3) = evecs_raw(absdof, md_i_sorted)* val_exp(ety_j)
            evecs_final_el(3,nd_p3) = evecs_raw(absdof, md_i_sorted)
         endif


      enddo


      !  Add the contribution of the element i_el to the mode component
      mode_poln_fracs(1:3, md_i) = mode_poln_fracs(1:3, md_i) + mode_comp(1:3)


      ! Copy this element into the main evecs_finalution array
      evecs_final(:,:,md_i,i_el) =  evecs_final_el

   enddo  ! end of current element




   !  Total energy and normalization
   z_tmp2 = mode_poln_fracs(1,md_i) + mode_poln_fracs(2,md_i) + mode_poln_fracs(3,md_i)
   if (abs(z_tmp2) < 1.0d-20) then ! 11/12/2024, trying to allow thin triangle element

      write(emsg,*) "The total energy for mode ", md_i, "is too small : ", z_tmp2
      call nberr%set(NBERR_BAD_ELT_ENERGY, emsg)
      return
   endif

   mode_poln_fracs(:,md_i) = mode_poln_fracs(:,md_i) / z_tmp2

   !  Check if the eigenvector is nonzero
   if (abs(z_evecs_final_max) < 1.0d-20) then ! 11/12/2024, trying to allow thin triangle element
      write(emsg,*) "The largest node value for mode ", md_i, "is too small : ", z_evecs_final_max
      call nberr%set(NBERR_BAD_ELT_ENERGY, emsg)
      return
   endif

   !  Normalization for this mode so that the maximum field component has magnitude 1
   evecs_final(:,:,md_i,:) = evecs_final(:,:,md_i,:)/z_evecs_final_max
   evecs_raw(1:cscmat%n_dof,md_i_sorted) = evecs_raw(1:cscmat%n_dof,md_i_sorted)/z_evecs_final_max

end subroutine









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
   use numbatmod

   integer(8), intent(in) :: n_modes
   complex(8), intent(in) :: shift_ksqr
   complex(8) :: v_evals_beta(n_modes)
   integer(8) :: v_eig_index(n_modes)

   integer(8) i

   complex(8) z_beta
   type(NBError) nberr

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
   call find_eigvals_order (n_modes, v_evals_beta, v_eig_index, 0_8, nberr)

   ! Apply the reordering
   call zvec_reorder_by_index(v_evals_beta, v_eig_index, n_modes)

end subroutine


! subroutine make_pbc_phase_shifts(mesh, entities, pbcs, i_el, bloch_vec, val_exp)

!    use numbatmod
!    use class_Mesh
!    use class_PeriodicBCs

!    type(MeshEM) :: mesh
!    type(MeshEntities) :: entities
!    type(PeriodicBCs) :: pbcs
!    double precision bloch_vec(2)

!    complex(8) val_exp(N_ENTITY_PER_EL)

!    integer(8) i_el, nd_i, mesh_pt, j1, ety_tag, j, k

!    integer(8) el_nds_i(P2_NODES_PER_EL)



!    complex(8) r_tmp1
!    double precision delta_xy_ref(2)
!    double precision ddot

!    val_exp = D_ONE

!    do nd_i=1,P2_NODES_PER_EL
!       mesh_pt = mesh%m_elnd_to_mshpt(nd_i,i_el)
!       k = pbcs%iperiod_N(mesh_pt)
!       if (k /= 0) mesh_pt=k
!       el_nds_i(nd_i) = mesh_pt
!    enddo

!    !  val_exp: Bloch mod ephase factor between the origin point and destination point
!    !  For a pair of periodic points, one is chosen as origin and the other is the destination
!    do j=1,N_ENTITY_PER_EL
!       ety_tag = entities%v_tags(j,i_el)
!       j1 = pbcs%iperiod_N_E_F(ety_tag)
!       if (j1 /= 0) then
!          !do k=1,dim_32
!          !  delta_xy_ref(k) = entities.v_mshpt_xy(k,ety_tag) - entities.v_mshpt_xy(k,j1)
!          !enddo
!          delta_xy_ref = entities%v_xy(:,ety_tag) - entities%v_xy(:,j1)
!          r_tmp1 = ddot(2, bloch_vec, 1, delta_xy_ref, 1)
!          val_exp(j) = exp(C_IM_ONE * r_tmp1)
!       endif
!    enddo

! end subroutine

