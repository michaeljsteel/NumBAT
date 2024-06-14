
#include "numbat_decl.h"

module nbinterfaces


    interface
        subroutine prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, &
                int_max, cmplx_max, real_max, &
                awk, bwk, cwk, dwk, iindex, overlap_L, &
                errco, emsg)

            use numbatmod

            integer*8 is_em, n_msh_el, n_msh_pts, n_modes
            integer*8 int_max, cmplx_max, real_max

            integer*8, dimension(:), allocatable :: awk
            complex*16, dimension(:), allocatable :: bwk
            double precision, dimension(:), allocatable :: cwk
            double precision, dimension(:,:), allocatable :: dwk
            integer*8, dimension(:), allocatable :: iindex
            complex*16, dimension(:,:), allocatable :: overlap_L

            integer*8 errco
            character(len=EMSG_LENGTH) emsg

        end subroutine prepare_workspaces

    end interface


end module nbinterfaces
