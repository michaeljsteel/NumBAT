
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


end module nbinterfaces
