
module nbinterfaces

    interface
        subroutine make_csr_arrays(mesh_raw, entities,  neq,&
            m_eqs, nonz, v_row_ind, v_col_ptr, debug, errco, emsg)

            use numbatmod
            use alloc
            use class_MeshRaw

            type(MeshRaw) :: mesh_raw
            type(MeshEntities) :: entities

            integer(8) neq, nonz
            integer(8) m_eqs(3, entities%n_ddl)
            !integer(8) m_eqs(:,:)


            integer(8), dimension(:), allocatable, intent(inout) :: v_row_ind
            integer(8), dimension(:), allocatable, intent(inout) :: v_col_ptr

            integer(8) debug

            integer(8) errco
            character(len=EMSG_LENGTH) emsg

        end subroutine

    end interface


end module nbinterfaces

module nbinterfacesb
    use numbatmod
    implicit none


    interface

        subroutine csr_length (nel, n_ddl, neq,  &
                table_N_E_F, ineq, &
                ext_v_row_ind, ext_v_col_ptr, &
                nonz_max, nonz, max_row_len, debug, errco, emsg)


            integer(8) nel, n_ddl, neq
            integer(8) table_N_E_F(14,nel)
            integer(8) ineq(3,n_ddl)

            integer(8), dimension(:), allocatable, intent(inout) :: ext_v_row_ind
            integer(8), dimension(:) :: ext_v_col_ptr(neq+1)

            integer(8) nonz_max, nonz, max_row_len

            integer(8) debug

            integer(8) errco
            !character(len=EMSG_LENGTH) emsg
            character(len=*) emsg
        end subroutine

    end interface


end module nbinterfacesb

