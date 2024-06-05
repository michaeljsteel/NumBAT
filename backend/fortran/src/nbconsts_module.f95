module nbconsts
      implicit none

      double precision, parameter :: D_ONE = 1.0d0
      double precision, parameter :: D_ZERO = 0.0d0
      complex*16, parameter :: C_IM_ONE = cmplx(0.0d0, 1.0d0, 8)

      integer*8, parameter :: nnodes_0 = 6
      integer*8, parameter :: nddl_t = 4


      integer*8, parameter :: BCS_CLOSED = 0
      integer*8, parameter :: BCS_OPEN = 1
      integer*8, parameter :: BCS_PERIODIC = 2

end module
