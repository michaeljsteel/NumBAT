
module nbinterfaces


   interface
      subroutine prepare_workspaces(n_msh_pts, n_msh_el, n_modes, int_max, cmplx_max, real_max, &
         awk, bwk, cwk, overlap_L, iindex, errco, emsg)

      integer*8 int_max, cmplx_max, real_max
      integer*8 n_msh_el, n_msh_pts, n_modes

      integer*8, dimension(:), allocatable :: awk
      complex*16, dimension(:), allocatable :: bwk
      double precision, dimension(:), allocatable :: cwk
      integer*8, dimension(:), allocatable :: iindex
      complex*16, dimension(:,:), allocatable :: overlap_L
      integer*8 errco
      character*2048 emsg
      end subroutine prepare_workspaces

   end interface


end module nbinterfaces
