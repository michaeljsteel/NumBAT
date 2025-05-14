#include "numbat_decl.h"

!  Difference from array_evecs_final_AC.f is that the u_z fi_eld is multiplied by i
!  which gives you the correct physical displacement fi_eld.

!  evecs_raw(*,i) : contains the imaginary and real parts of the solution for points such that cscmat%m_global_dofs(i) != 0
!  evecs_final(i) : contains solution for all points
!  The dimension of the geometric domain is : dim_32 = 2
!  The dimension of the vector fi_eld is : dim2 = 3

! Mode polarisation fractions are currently turned off

subroutine construct_solution_fields_ac (shift_omsq, n_modes, mesh, cscmat,  &
   v_evals_nu, evecs_raw, evecs_final, mode_poln_fracs, nberr)


   use numbatmod
   use alloc

   use class_Mesh
   use class_SparseCSC_AC

   type(MeshAC) mesh
   type(SparseCSC_AC) cscmat
   type(NBError) nberr


   double precision shift_omsq
   integer(8) n_modes

   !  TODO: n_core seems to be never initialised. Is that code ever called?
   !integer(8) n_core(2)

   integer(8), dimension(:), allocatable :: v_eig_index

   complex(8) evecs_raw(cscmat%n_dof,n_modes)

   !  evecs_final(3, 1..P2_NODES_PER_EL,n_modes, mesh%n_msh_elts)
   ! contains the values of the 3 components at P2 interpolation nodes
   complex(8) evecs_final(3,P2_NODES_PER_EL,n_modes,mesh%n_msh_elts)
   complex(8) v_evals_nu(n_modes) !, v_tmp(n_modes)
   complex(8) mode_poln_fracs(4,n_modes)



   double precision mode_comp(4)
   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2,P2_NODES_PER_EL)
   complex(8) evecs_final_el(3,P2_NODES_PER_EL)


   integer(8) nd_i, typ_e, debug
   integer(8) i_el, md_i, md_i2, msh_pt_i, dof, nd_xyz
   complex(8) z_tmp2, z_evecs_final_max, z_beta

   character(len=EMSG_LENGTH) :: emsg

   debug = 0


   ! rescale and sort eigenvalues

   call integer_alloc_1d(v_eig_index, n_modes, 'v_eig_index', nberr); RET_ON_NBERR(nberr)

   do md_i=1,n_modes
      z_beta = sqrt(1.0d0/v_evals_nu(md_i) + shift_omsq) / (2.0d0 * D_PI)
      !!  Frequency (z_beta) should always be positive.
      v_evals_nu(md_i) = abs(z_beta)
   enddo

   ! order with smallest eigenvalue first (order=1)
   call find_eigvals_order (n_modes, v_evals_nu, v_eig_index, 1_8, nberr)
   call zvec_reorder_by_index(v_evals_nu, v_eig_index, n_modes)


   ! reorder and normalise the eigenvectors

   evecs_final= C_ZERO
   z_tmp2=D_ZERO
   mode_poln_fracs = D_ZERO

   do md_i=1,n_modes
      md_i2 = v_eig_index(md_i) !  index of the mode in eigenvalue sorted sequence


      z_evecs_final_max = 0.0d0 !  value and loc of max field modulus

      do i_el=1,mesh%n_msh_elts            ! for each elt
         typ_e = mesh%v_elt_material(i_el)

         mode_comp = D_ZERO

         call mesh%find_nodes_for_elt(i_el, el_nds_i, el_nds_xy)

         do nd_i=1,P2_NODES_PER_EL           ! for each of 6 nodes of that elt
            msh_pt_i = mesh%m_elnd_to_mshpt(nd_i,i_el)  ! map to the actual mesh pt

            do nd_xyz=1,3                    ! for each xyz component of this mesh pt
               dof = cscmat%m_global_dofs(nd_xyz, msh_pt_i)  ! the actual fem dof
               if (dof .gt. 0) then
                  evecs_final_el(nd_xyz, nd_i) = evecs_raw(dof, md_i2)
               else
                  evecs_final_el(nd_xyz, nd_i) = 0
               endif
            enddo

            !  The z-component must be multiplied by ii in order to get the un-normalised z-component
            nd_xyz=3
            evecs_final_el(nd_xyz,nd_i) = C_IM_ONE* evecs_final_el(nd_xyz,nd_i)

            do nd_xyz=1,3
               z_tmp2 = evecs_final_el(nd_xyz, nd_i)
               evecs_final(nd_xyz, nd_i, md_i, i_el) = z_tmp2
               if (abs(z_evecs_final_max) .lt. abs(z_tmp2)) then
                  z_evecs_final_max = z_tmp2
                  !  We want to normalise such the the z-component is purely imaginary complex number
                  if (nd_xyz == 3) z_evecs_final_max = - C_IM_ONE* z_evecs_final_max
               endif
            enddo

            !  Contribution of the element i_el to the mode component
            ! mode_comp(1:3) = mode_comp(1:3) + abs(evecs_final_el(:,nd_i))**2

         enddo


         !  Add the contribution of the element i_el to the mode component
         ! mode_poln_fracs(1:3,md_i) = mode_poln_fracs(1:3,md_i) + mode_comp(1:3)/P2_NODES_PER_EL

         ! !  TODO: THIS CODE SEEMS NEVER TO BE CALLED
         ! !  THESE BAD VALUES TO TRY TO MAKE IT FAIL IF IT DOES
         ! n_core(1)=-1000
         ! n_core(2)=-1000

         ! if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
         !    write(*,*) "Warning: unitialised n_core in array_evecs_final_AC.f"
         !    z_tmp2 = mode_comp(1) + mode_comp(2)&
         !    &+ mode_comp(3)
         !    mode_poln_fracs(4,md_i) = mode_poln_fracs(4,md_i) + z_tmp2
         ! endif

      enddo

      ! !  Total energy and normalization
      ! z_tmp2 = mode_poln_fracs(1,md_i) + mode_poln_fracs(2,md_i) + mode_poln_fracs(3,md_i)
      ! if (abs(z_tmp2) .lt. 1.0d-10) then
      !    write(emsg,*) "array_evecs_final: the total energy ", "is too small : ", z_tmp2, &
      !     "array_evecs_final: md_i md_i2 = ", md_i, md_i2
      !    call nberr%set(NBERR_BAD_ELASTIC_ENERGY, emsg)
      !    return
      ! endif

      !    mode_poln_fracs(:,md_i) = mode_poln_fracs(:,md_i) / z_tmp2


      !  Check if the eigenvector is nonzero
      if (abs(z_evecs_final_max) .lt. 1.0d-10) then
         z_evecs_final_max = z_tmp2
         write(emsg,*) "array_evecs_final: z_evecs_final_max is too small" , z_evecs_final_max, &
            "at md_i, md_i2, n_modes = ", md_i, md_i2, n_modes
         call nberr%set(NBERR_BAD_ELASTIC_ENERGY, emsg)
         return
      endif


      ! !  Normalization so that the maximum field component is 1
      ! do i_el=1,mesh%n_msh_elts
      !    do nd_i=1,P2_NODES_PER_EL
      !       i1 = mesh%m_elnd_to_mshpt(nd_i,i_el)



      !    enddo
      ! enddo

      !evecs_raw(:,md_i2)  = evecs_raw(:,md_i2)/z_evecs_final_max

      evecs_final(:,:,md_i,:) = evecs_final(:,:,md_i,:) / z_evecs_final_max

   enddo


end
