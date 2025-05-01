
! module nbinterfaces

!     interface
!         subroutine make_csc_arraysb(mesh_raw, entities, cscmat, &
!             !v_row_ind, v_col_ptr,
!             debug, errco, emsg)

!             use numbatmod
!             use alloc
!             use class_MeshRawEM
!             use class_SparseCSC_EM

!             type(MeshRawEM) :: mesh_raw
!             type(MeshEntities) :: entities
!             type(SparseCSC_EM) :: cscmat

!             !integer(8), dimension(:), allocatable, intent(inout) :: v_row_ind
!             !integer(8), dimension(:), allocatable, intent(inout) :: v_col_ptr

!             integer(8) debug

!             integer(8) errco
!             character(len=EMSG_LENGTH) emsg

!         end subroutine

!     end interface


! end module nbinterfaces

module nbinterfacesb
    use numbatmod
    implicit none


    interface

        subroutine csr_length_b (nel, n_ddl, n_dof,  &
                table_N_E_F, in_dof, &
                ext_v_row_ind, ext_v_col_ptr, &
                nonz_max, nonz, max_row_len, errco, emsg)


            integer(8) nel, n_ddl, n_dof
            integer(8) table_N_E_F(14,nel)
            integer(8) in_dof(3,n_ddl)

            integer(8), dimension(:), allocatable, intent(inout) :: ext_v_row_ind
            integer(8), dimension(:) :: ext_v_col_ptr(n_dof+1)

            integer(8) nonz_max, nonz, max_row_len

            integer(8) errco
            !character(len=EMSG_LENGTH) emsg
            character(len=*) emsg
        end subroutine

    end interface


end module nbinterfacesb

