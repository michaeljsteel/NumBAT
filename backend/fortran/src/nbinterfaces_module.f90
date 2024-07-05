
module nbinterfaces

    use numbatmod
    implicit none

    interface
        subroutine prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, &
                int_max, cmplx_max, real_max, &
                a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L, &
                errco, emsg)


            integer, intent(in) :: is_em
            integer(8), intent(in) :: n_msh_el, n_msh_pts, n_modes
            integer(8), intent(out) :: int_max, cmplx_max, real_max

            integer(8), dimension(:), allocatable, intent(inout) :: a_iwork
            complex(8), dimension(:), allocatable, intent(inout) :: b_zwork
            double precision, dimension(:), allocatable, intent(inout) :: c_dwork
            double precision, dimension(:,:), allocatable, intent(inout) :: d_dwork
            integer(8), dimension(:), allocatable, intent(inout) :: iindex
            complex(8), dimension(:,:), allocatable, intent(inout) :: overlap_L

            integer, intent(out) :: errco
            character(len=*), intent(out) :: emsg

        end subroutine prepare_workspaces

    end interface

    interface

        subroutine csr_length (nel, n_ddl, neq,  &
                table_N_E_F, ineq, &
               ! col_ind, row_ptr, &
                ext_v_row_ind, ext_v_col_ptr, &
                nonz_max, nonz, max_row_len, int_max, debug, errco, emsg)


            integer(8) nel, n_ddl, neq
            integer(8) table_N_E_F(14,nel)
            integer(8) ineq(3,n_ddl)

            !integer(8) col_ind(*)
            !integer(8) row_ptr(neq+1)

            integer(8), dimension(:), allocatable, intent(inout) :: ext_v_row_ind
            integer(8), dimension(:) :: ext_v_col_ptr(neq+1)


            integer(8) nonz_max, nonz, max_row_len
            integer(8) int_max

            integer(8) debug

            integer errco
            !character(len=EMSG_LENGTH) emsg
            character(len=*) emsg
        end subroutine

    end interface


end module nbinterfaces
