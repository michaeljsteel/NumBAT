
#include "numbat_decl.h"

module numbatmod

   !  Use intel compiler to check passing conventions
#ifdef __INTEL_COMPILER
   use ifport
#endif

   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
      stdout=>output_unit, &
      stderr=>error_unit

   implicit none

   integer, parameter :: EMSG_LENGTH = 2048
   integer, parameter :: FNAME_LENGTH = 1024

   integer, parameter :: MAX_N_PTS = 250000
   integer, parameter :: MAX_N_ELTS = 120000

   integer, parameter :: MAX_LONG_ADJ = 2500000

   double precision, parameter :: D_PI = 3.141592653589793d0
   double precision, parameter :: D_ONE = 1.0d0
   double precision, parameter :: D_ZERO = 0.0d0

   complex(8), parameter :: C_ZERO = (0.0d0, 0.0d0)
   complex(8), parameter :: C_ONE = (1.0d0, 0.0d0)
   complex(8), parameter :: C_IM_ONE = (0.0d0, 1.0d0)


   double precision, parameter :: SI_C_SPEED = 299792458.0d0
   double precision, parameter :: SI_EPS_0 = 8.8541878188d-12
   double precision, parameter :: SI_MU_0 = 1.25663706127d-6




   integer(8), parameter :: nnodes_0 = 6
   integer(8), parameter :: nddl_t = 4

   integer(8), parameter :: nddl_0_em = 14
   integer(8), parameter :: nddl_0_ac = 6


   integer(8), parameter :: BCS_DIRICHLET = 0  !  (E-field: electric wall)
   integer(8), parameter :: BCS_NEUMANN = 1    !  (E-field: magnetic wall)
   integer(8), parameter :: BCS_PERIODIC = 2

   integer(8), parameter :: FEM_FORMULATION_E = 1
   integer(8), parameter :: FEM_FORMULATION_H = 2



   integer, parameter :: UMFPACK_CONTROL = 20
   integer, parameter :: UMFPACK_INFO = 90
   integer(8), parameter :: NBERR_BAD_PERMUTATION     = -52
   integer(8), parameter :: NBERR_BAD_ADJACENCY       = -53
   integer(8), parameter :: NBERR_BAD_NODE_SEPARATION = -54


   !  UMFPACK Solve codes
   !  Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the
   !  linear algebraic transpose (complex conjugate if A is complex), or the (')
   !  operator in MATLAB.  "at" refers to the array transpose, or the (.')
   !  operator in MATLAB.


   integer(8), parameter ::  UMFPACK_A      = 0     !  Ax=b
   integer(8), parameter ::  UMFPACK_At     = 1     !  A'x=b
   integer(8), parameter ::  UMFPACK_Aat    = 2     !  A.'x=b

   integer(8), parameter ::  UMFPACK_Pt_L    =3     !  P'Lx=b
   integer(8), parameter ::  UMFPACK_L       =4     !  Lx=b
   integer(8), parameter ::  UMFPACK_Lt_P    =5     !  L'Px=b
   integer(8), parameter ::  UMFPACK_Lat_P   =6     !  L.'Px=b
   integer(8), parameter ::  UMFPACK_Lt      =7     !  L'x=b
   integer(8), parameter ::  UMFPACK_Lat     =8     !  L.'x=b

   integer(8), parameter ::  UMFPACK_U_Qt    =9     !  UQ'x=b
   integer(8), parameter ::  UMFPACK_U       =10    !  Ux=b
   integer(8), parameter ::  UMFPACK_Q_Ut    =11    !  QU'x=b
   integer(8), parameter ::  UMFPACK_Q_Uat   =12    !  QU.'x=b
   integer(8), parameter ::  UMFPACK_Ut      =13    !  U'x=b
   integer(8), parameter ::  UMFPACK_Uat     =14    !  U.'x=b



contains

   subroutine log_me(fname, msg, newfile)
      character(len=*) :: fname, msg
      logical, optional :: newfile
      logical app

      integer(4) ui

      ui=8
      app = .false.
      if (present(newfile)) then
         app = .not. newfile
      endif

      if (app) then
         open(ui, FILE=fname, ACTION='WRITE', ACCESS='APPEND')
      else
         open(ui, FILE=fname, ACTION='WRITE')
      endif

      write(ui, *) msg

      close(ui)

   end subroutine

   subroutine flush_me(ui, msg)
      ! Declare the interface for POSIX fsync function
      interface
         function fsync (fd) bind(c,name="fsync")
            use iso_c_binding, only: c_int
            integer(c_int), value :: fd
            integer(c_int) :: fsync
         end function fsync
      end interface

      integer(8) :: ui, ret
      character(len=*) :: msg

      write(ui, '(A,A)') '>>>> ', msg
      flush(ui)

#ifndef __INTEL_COMPILER
            ret = fsync(fnum(ui))   ! actually sync to fs on GCC
#endif

   end subroutine

   integer function nb_system(cmd)
      implicit none

      integer errco
      character(len=*), intent(in) :: cmd

#ifdef __INTEL_COMPILER
      errco = system(cmd)
#else
      call system(cmd, errco)
#endif

      nb_system = errco

   end function

   logical function almost_equal(a, b) result(res)
      implicit none

      real(8) :: a, b
      real(8), parameter :: tol=1.d-12

      res = abs(a-b) < tol

   end function

   subroutine assert_or_die(pred, msg, ec)  !  TODO: this is cheat rather than go back to python for reporting
      logical :: pred
      character(len=*) :: msg
      integer :: ec

      if (pred) return

      write(*,*) msg
      call exit(ec)

   end subroutine

   subroutine assert_no_larger_than(val, limit, location, msg, failco, errco, emsg)

      implicit none

      integer errco
      character(len=EMSG_LENGTH) emsg
      character location*(*), msg*(*)
      integer val, limit, failco

      if (val .ge. limit) then
         write(emsg,*) 'Failed limit check at ', location, '.  ', &
            'Expected ', msg, ',  but found values', val, limit
         errco = failco
      endif

      return
   end subroutine

   function int_2_str(val, fmt) result(str)
      integer(8), intent(in) :: val
      character(len=*), intent(in), optional :: fmt

      character(len=:), allocatable :: str

      integer, parameter :: buflen = 512
      character(len=buflen) :: buffer

      character(len=buflen) :: d_fmt = '(i0)'

      if (present(fmt)) then
         d_fmt = fmt
      endif

      write(buffer, d_fmt) val

      str = trim(buffer)
   end function int_2_str

   !  TODO: fix with just one call
   function int4_2_str(val, fmt) result(str)
      integer(4), intent(in) :: val
      character(len=*), intent(in), optional :: fmt

      character(len=:), allocatable :: str

      integer, parameter :: buflen = 512
      character(len=buflen) :: buffer

      character(len=buflen) :: d_fmt = '(i0)'

      if (present(fmt)) then
         d_fmt = fmt
      endif

      write(buffer, d_fmt) val

      str = trim(buffer)
   end function int4_2_str

   function double_2_str(val, fmt) result(str)
      double precision, intent(in) :: val
      character(len=*), intent(in), optional ::fmt

      character(len=:), allocatable :: str
      integer, parameter :: buflen = 512
      character(len=buflen) :: buffer

      character(len=buflen) :: d_fmt = '(e)'

      if (present(fmt)) then
         d_fmt = fmt
      endif

      write(buffer, d_fmt) val

      str = trim(buffer)
   end function double_2_str


   subroutine get_clocks(systime, cputime)
      !  Returns system (wall time) in seconds, and cpu time in seconds
      !  nanosec may be microsec on some systems

      integer(8) isystime
      double precision systime, cputime
      double precision nanosec

      parameter (nanosec=1.d-9)

      call system_clock(isystime)
      call cpu_time(cputime)

      systime = nanosec*isystime

   end subroutine get_clocks



   !TODO: move somewhere leass general
   logical function log_is_curved_elem_tri (nnodes, xel) result(is_curved)

      implicit none
      integer(8) nnodes, info_curved
      double precision xel(2,nnodes)

      double precision loctmp

      call is_curved_elem_tri_impl (nnodes, xel, info_curved, loctmp)

      is_curved = (info_curved .ne. 0)   !  painful, eventually get rid of the int form of is_curved_elem_tri


   end function

end module numbatmod
