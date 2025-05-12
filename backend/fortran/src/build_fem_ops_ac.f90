
#include "numbat_decl.h"
subroutine build_fem_ops_ac (shift_omsq, q_ac, rho, c_tensor, &
   mesh_raw, cscmat, symmetry_flag, nberr)


   use numbatmod
   use alloc
   use class_Mesh
   use class_SparseCSC_AC
   use class_BasisFunctions

   type(MeshAC) mesh_raw
   type(SparseCSC_AC) cscmat
   type(NBError) nberr


   integer(8) symmetry_flag
   complex(8) shift_omsq, q_ac

   complex(8) rho(mesh_raw%n_elt_mats), c_tensor(6,6,mesh_raw%n_elt_mats)


   integer(8), dimension(:), allocatable :: i_work

   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2, P2_NODES_PER_EL)

   integer(8) typ_e
   integer(8) i, k, el_i
   integer(8) nd_i, msh_pt_i, dof_i, nd_dof_i
   integer(8) nd_j, msh_pt_j, dof_j, nd_dof_j
   integer(8) col_start, col_end

   complex(8) mat_K(18,18), mat_M(18,18)
   complex(8) c_tensor_el(6,6), rho_el
   complex(8) K_elt, M_elt

   type(BasisFunctions) basfuncs

   character(len=EMSG_LENGTH) :: emsg


   call integer_alloc_1d(i_work, 3*mesh_raw%n_msh_pts, 'i_work', nberr); RET_ON_NBERR(nberr)


   cscmat%mOp_stiff = C_ZERO
   cscmat%mOp_mass = C_ZERO

   do el_i=1,mesh_raw%n_msh_elts              ! For each elt
      typ_e = mesh_raw%v_elt_material(el_i)    ! Find the material and its local material properties

      rho_el = rho(typ_e)
      c_tensor_el = c_tensor(:,:,typ_e)

      call mesh_raw%find_nodes_for_elt(el_i, el_nds_i, el_nds_xy)

      call basfuncs%set_affine_for_elt(el_nds_xy, nberr)
      RET_ON_NBERR(nberr)
      call basfuncs%get_triint_set_p2_p2()


      !  If c_tensor has regular symmetries use more efficient formulation
      !write(*,*) 'symflag', symmetry_flag
      if (symmetry_flag .eq. 1) then
         call mat_el_v2 (el_nds_xy, q_ac, c_tensor_el, rho_el, mat_K, mat_M)
      else
         call make_elt_femops_ac_no_sym (basfuncs, q_ac, c_tensor_el, rho_el, mat_K, mat_M)
      endif

      do nd_j=1,P2_NODES_PER_EL  ! iterating columns
         msh_pt_j = mesh_raw%m_elnd_to_mshpt(nd_j,el_i)

         do nd_dof_j=1,3
            dof_j = cscmat%m_eqs(nd_dof_j,msh_pt_j)

            if (dof_j .gt. 0) then
               col_start = cscmat%v_col_ptr(dof_j)
               col_end = cscmat%v_col_ptr(dof_j+1) - 1

               !  unpack row into i_work
               do i=col_start,col_end
                  i_work(cscmat%v_row_ind(i)) = i
               enddo

               do nd_i=1,P2_NODES_PER_EL   ! iterating rows
                  msh_pt_i = mesh_raw%m_elnd_to_mshpt(nd_i,el_i)

                  do nd_dof_i=1,3
                     dof_i = cscmat%m_eqs(nd_dof_i,msh_pt_i)

                     if (dof_i .gt. 0) then
                        K_elt = mat_K(3*(nd_i-1) + nd_dof_i, 3*(nd_j-1) + nd_dof_j)
                        M_elt = mat_M(3*(nd_i-1) + nd_dof_i, 3*(nd_j-1) + nd_dof_j)
                        K_elt = K_elt - shift_omsq*M_elt
                        k = i_work(dof_i)
                        if (k .gt. 0 .and. k .le. cscmat%n_nonz) then

                           cscmat%mOp_stiff(k) = cscmat%mOp_stiff(k) + K_elt
                           cscmat%mOp_mass(k) = cscmat%mOp_mass(k) + M_elt

                        else
                           write(emsg,*) "build_fem_ops_ac: problem with  row_ind !!", k, cscmat%n_nonz
                           call nberr%set(NBERR_BAD_ASSEMBLY_AC, emsg);
                        endif
                     endif

                  enddo
               enddo
            endif
         enddo
      enddo
   enddo

end
