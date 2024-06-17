
module nbinterfaces


   interface
      subroutine prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, &
         int_max, cmplx_max, real_max, &
         a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L, &
         errco, emsg)

         use numbatmod

         integer :: is_em, n_modes
         integer*8 :: n_msh_el, n_msh_pts
         integer*8 :: int_max, cmplx_max, real_max


         integer*8, dimension(:), allocatable :: a_iwork
         complex*16, dimension(:), allocatable :: b_zwork
         double precision, dimension(:), allocatable :: c_dwork
         double precision, dimension(:,:), allocatable :: d_dwork
         integer*8, dimension(:), allocatable :: iindex
         complex*16, dimension(:,:), allocatable :: overlap_L

         integer :: errco
         character(len=EMSG_LENGTH) emsg

      end subroutine prepare_workspaces

   end interface


end module nbinterfaces
