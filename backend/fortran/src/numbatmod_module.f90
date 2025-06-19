
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

#include "nbversion_incl.h"

   integer(8), parameter :: P1_NODES_PER_EL = 3
   integer(8), parameter :: P2_NODES_PER_EL = 6
   integer(8), parameter :: P3_NODES_PER_EL = 10


   ! elastic fem properties
   integer(8), parameter :: N_DOF_PER_NODE_AC = 3  ! xyz components for each node
   !integer(8), parameter :: N_DOF_PER_NODE_AC_PIEZO = 6  ! U_xyz and D_xyz components for each node


   ! sequencing of the P3 nodes

   ! indices into MeshEntities%v_xy
   !   ** these elements start at these values +1 **
   integer(8), parameter :: ETY_TAG_OFFSET_FACE = 0         ! uses slot 1
   integer(8), parameter :: ETY_TAG_OFFSET_P2_EDGES = 1     ! uses slots 2,3,4
   integer(8), parameter :: ETY_TAG_OFFSET_P3_VERTICES = 4  ! uses slots 5,6,7
   integer(8), parameter :: ETY_TAG_OFFSET_P3_NODES = 7     ! uses slots 8..13
   integer(8), parameter :: ETY_TAG_OFFSET_P3_INTERIOR = 13 ! uses slots 8..14

   ! indices into MeshEntities%v_ety_props
   integer(8), parameter :: ETY_PROP_PHYSTYPE = 1
   integer(8), parameter :: ETY_PROP_DIMENSION = 2




   integer(8), parameter :: P3_VERT_1 = 1
   integer(8), parameter :: P3_VERT_2 = 2
   integer(8), parameter :: P3_VERT_3 = 3
   integer(8), parameter :: P3_EDGE_LO = 4
   integer(8), parameter :: P3_EDGE_HI = 9
   integer(8), parameter :: P3_INTERIOR = 10


   integer(8), parameter :: N_ETY_TRANSVERSE = 4   ! Number of transverse dof
   integer(8), parameter :: N_ENTITY_PER_EL = 14   ! 1 Face + 6 P2 nodes + 7 additional P3 nodes
   integer(8), parameter :: N_DOF_PER_EL = 13      ! 6 P2 nodes + 7 additional P3 nodes


!   integer(8), parameter :: NDDL_0_AC = 6

   integer(8),  parameter :: EMSG_LENGTH = 2048
   integer(8),  parameter :: FNAME_LENGTH = 1024

   integer(8), parameter :: MAX_N_PTS = 250000
   integer(8), parameter :: MAX_N_ELTS = 120000

   integer(8), parameter :: MAX_LONG_ADJ = 2500000

   double precision, parameter :: D_PI = 3.141592653589793d0
   double precision, parameter :: D_ONE = 1.0d0
   double precision, parameter :: D_ZERO = 0.0d0

   complex(8), parameter :: C_ZERO = (0.0d0, 0.0d0)
   complex(8), parameter :: C_ONE = (1.0d0, 0.0d0)
   complex(8), parameter :: C_IM_ONE = (0.0d0, 1.0d0)


   double precision, parameter :: SI_C_SPEED = 299792458.0d0
   double precision, parameter :: SI_EPS_0 = 8.8541878188d-12
   double precision, parameter :: SI_MU_0 = 1.25663706127d-6




   integer(8), parameter :: BCS_DIRICHLET = 0  !  (E-field: electric wall)
   integer(8), parameter :: BCS_NEUMANN = 1    !  (E-field: magnetic wall)
   integer(8), parameter :: BCS_PERIODIC = 2

   integer(8), parameter :: FEM_FORMULATION_E = 1
   integer(8), parameter :: FEM_FORMULATION_H = 2



   integer(8),  parameter :: UMFPACK_CONTROL = 20
   integer(8),  parameter :: UMFPACK_INFO = 90
   integer(8), parameter :: NBERR_BAD_PERMUTATION     = -52
   integer(8), parameter :: NBERR_BAD_ADJACENCY       = -53
   integer(8), parameter :: NBERR_BAD_NODE_SEPARATION = -54
   integer(8), parameter :: NBERR_BAD_JACOBIAN        = -55
   integer(8), parameter :: NBERR_BAD_DETERMINANT     = -56
   integer(8), parameter :: NBERR_BAD_ASSEMBLY        = -57
   integer(8), parameter :: NBERR_BAD_MESH_EDGES      = -58
   integer(8), parameter :: NBERR_BAD_MESH_VERTICES   = -59
   integer(8), parameter :: NBERR_BAD_BOUNDARY_CONDITION   = -60
   integer(8), parameter :: NBERR_BAD_QUAD_INT  = -61
   integer(8), parameter :: NBERR_BAD_MB_EDGES  = -62
   integer(8), parameter :: NBERR_BAD_ELT_ENERGY  = -63
   integer(8), parameter :: NBERR_MESH_TOO_LARGE  = -64
   integer(8), parameter :: NBERR_SORT_STACK_TOO_SMALL  = -65
   integer(8), parameter :: NBERR_UNKNOWN_SORT_ORDER = -66
   integer(8), parameter :: NBERR_BAD_ELASTIC_ENERGY = -67
   integer(8), parameter :: NBERR_BAD_ASSEMBLY_AC = -68
   integer(8), parameter :: NBERR_BAD_ZNAUPD = -106
   integer(8), parameter :: NBERR_BAD_UMF4ZSOL = -107
   integer(8), parameter :: NBERR_BAD_UMF4ZSYM = -108
   integer(8), parameter :: NBERR_BAD_UMF4ZNUM = -109
   integer(8), parameter :: NBERR_BAD_SPARSE_FINAL = -110




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



   type, public :: NBError

      integer(8) errco
      character(len=EMSG_LENGTH) :: emsg
   contains

      procedure :: reset => NBError_reset
      procedure :: ok => NBError_ok
      procedure :: set => NBError_set
      procedure :: to_py => NBError_to_py

   end type NBError


contains

   function NBError_ok (this) result(is_ok)
      class(NBError) this
      logical is_ok
      is_ok = (this%errco .eq. 0)
   end function

   subroutine NBError_reset (this)
      class(NBError) this
      this%errco = 0
      this%emsg = ""
   end subroutine


   subroutine NBError_set (this, ec, msg)
      class(NBError) this
      integer(8) ec
      character(len=*) :: msg
      this%errco = ec
      write(this%emsg, '(A)') msg

   end subroutine

   subroutine NBError_to_py(this, ec, msg)
      class(NBError) this
      integer(8), intent(out) :: ec
      character(len=EMSG_LENGTH), intent(out) :: msg

      ec = this%errco
      write(msg, '(A)') this%emsg

   end subroutine



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

      character(len=*) :: msg
      integer(8) :: ui

      ! #ifdef __GNUC__
      !       integer(8) :: ret
      ! #endif

      write(ui, '(A,A)') '>>>> ', msg
      flush(ui)

      ! #ifdef __GNUC__
      !             ret = fsync(fnum(ui))   ! actually sync to fs on GCC
      ! #endif

   end subroutine

   integer(8) function nb_system(cmd)
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
      integer(8) :: ec

      if (pred) return

      write(*,*) msg
      call exit(int(ec,4))

   end subroutine

   subroutine assert_no_larger_than(val, limit, location, msg, failco, nberr)

      type(NBError) nberr

      character(len=EMSG_LENGTH) emsg
      character location*(*), msg*(*)
      integer(8) val, limit
      integer(8) failco

      if (val .ge. limit) then
         write(emsg,*) 'Failed limit check at ', location, '.  ', &
            'Expected ', msg, ',  but found values', val, limit
            call nberr%set(failco, emsg)
      endif

   end subroutine

   function int_2_str(val, fmt) result(str)
      integer(8), intent(in) :: val
      character(len=*), intent(in), optional :: fmt

      character(len=:), allocatable :: str

      integer(8),  parameter :: buflen = 512
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

      integer(8),  parameter :: buflen = 512
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
      integer(8),  parameter :: buflen = 512
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

   complex(8) function v2_dot_v3(dveca, cvecb)
      double precision dveca(2)
      complex (8) cvecb(3)
      v2_dot_v3 = dveca(1) * cvecb(1) +  dveca(2) * cvecb(2)
   end function

   complex(8) function cv3_dot_cv3(cveca, cvecb)
      complex (8) cveca(3), cvecb(3)
      cv3_dot_cv3 = cveca(1) * cvecb(1) + cveca(2) * cvecb(2) + cveca(3) * cvecb(3)
   end function



   subroutine v2_cross_v3(dveca, cvecb, vres)
      complex(8) vres(3)
      double precision dveca(2)
      complex (8) cvecb(3)


      vres(1) = dveca(2) * cvecb(3)   ! fx = gy hz - gz hy = gy hz
      vres(2) = -dveca(1) * cvecb(3)  ! fy = gz hx - gx hz = -gx hz
      vres(3) = dveca(1) * cvecb(2) - dveca(2) * cvecb(1)  ! fz = gx hy - gy hx
   end subroutine




end module numbatmod
