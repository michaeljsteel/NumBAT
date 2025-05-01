
#include "numbat_decl.h"
subroutine assembly_ac (i_base, shift_omsq, q_ac, rho, c_tensor, &
   mesh_raw, cscmat, symmetry_flag, nberr)


   use numbatmod
   use alloc
   use class_MeshRawEM
   use class_SparseCSC_AC

   type(MeshRawAC) mesh_raw
   type(SparseCSC_AC) cscmat
   type(NBError) nberr


   integer(8) i_base

   integer(8) symmetry_flag
   complex(8) shift_omsq, q_ac

   complex(8) rho(mesh_raw%n_elt_mats), c_tensor(6,6,mesh_raw%n_elt_mats)


   integer(8), dimension(:), allocatable :: i_work

   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2, P2_NODES_PER_EL)

   integer(8) i_base2,  typ_e
   integer(8) i, k, i_el
   integer(8) ety_j, msh_pt_j, eqn_j, dof_j
   integer(8) ety_i, msh_pt_i, eqn_i, dof_i
   integer(8) col_start, col_end

   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el
   complex(8) z_tmp1, z_tmp2

!  The CSC indexing, i.e., cscmat%v_col_ptr, is 1-based
!  But valpr.f may have changed the CSC indexing to 0-based indexing)

   call integer_nalloc_1d(i_work, 3*mesh_raw%n_msh_pts, 'i_work', nberr); RET_ON_NBERR(nberr)

   if (i_base .eq. 0) then
      i_base2 = 1
   else
      i_base2 = 0
   endif

   cscmat%mOp_stiff = C_ZERO
   cscmat%mOp_mass = C_ZERO

   do i_el=1,mesh_raw%n_msh_el              ! For each elt

      typ_e = mesh_raw%el_material(i_el)    ! Find the material and its local material properties

      rho_el = rho(typ_e)
      c_tensor_el = c_tensor(:,:,typ_e)


      call mesh_raw%find_nodes_for_elt(i_el, el_nds_i, el_nds_xy)

      !  If c_tensor has regular symmetries use more efficient formulation
      if (symmetry_flag .eq. 1) then
         call mat_el_v2 (el_nds_xy,q_ac,c_tensor_el,rho_el,mat_K,mat_M)
      else
         call mat_el_v3 (el_nds_xy,q_ac,c_tensor_el,rho_el,mat_K,mat_M)
      endif

      do ety_j=1,P2_NODES_PER_EL
         msh_pt_j = mesh_raw%elnd_to_mshpt(ety_j,i_el)

         do dof_j=1,3
            eqn_j = cscmat%m_eqs(dof_j,msh_pt_j)

            if (eqn_j .gt. 0) then
               col_start = cscmat%v_col_ptr(eqn_j) + i_base2
               col_end = cscmat%v_col_ptr(eqn_j+1) - 1 + i_base2

               !  unpack row into i_work
               do i=col_start,col_end
                  i_work(cscmat%v_row_ind(i) + i_base2) = i
               enddo

               do ety_i=1,P2_NODES_PER_EL
                  msh_pt_i = mesh_raw%elnd_to_mshpt(ety_i,i_el)

                  do dof_i=1,3
                     eqn_i = cscmat%m_eqs(dof_i,msh_pt_i)

                     if (eqn_i .gt. 0) then
                        z_tmp1 = mat_K(3*(ety_i-1) + dof_i, 3*(ety_j-1) + dof_j)
                        z_tmp2 = mat_M(3*(ety_i-1) + dof_i, 3*(ety_j-1) + dof_j)
                        z_tmp1 = z_tmp1 - shift_omsq*z_tmp2
                        k = i_work(eqn_i)
                        if (k .gt. 0 .and. k .le. cscmat%n_nonz) then

                           cscmat%mOp_stiff(k) = cscmat%mOp_stiff(k) + z_tmp1
                           cscmat%mOp_mass(k) = cscmat%mOp_mass(k) + z_tmp2

                        else
                           write(*,*) "asmbly_AC: problem with cscmat%v_row_ind !!"
                           write(*,*) "asmbly_AC: k, cscmat%n_nonz = ", k, cscmat%n_nonz
                           write(*,*) "asmbly_AC: Aborting..."
                           stop
                        endif
                     endif

                  enddo
               enddo
            endif
         enddo
      enddo
   enddo

end
