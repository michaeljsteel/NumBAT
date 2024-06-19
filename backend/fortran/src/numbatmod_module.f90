
#include "numbat_decl.h"

module numbatmod

    ! Use intel compiler to check passing conventions
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
   !complex(8), parameter :: C_IM_ONE = cmplx(0.0d0, 1.0d0, 8)
   complex(8), parameter :: C_IM_ONE = (0.0d0, 1.0d0)


   double precision, parameter :: SI_C_SPEED = 299792458.0d0
   double precision, parameter :: SI_EPS_0 = 8.8541878188d-12
   double precision, parameter :: SI_MU_0 = 1.25663706127d-6


   integer(8), parameter :: nodes_per_el = 6

   integer(8), parameter :: nnodes_0 = 6
   integer(8), parameter :: nddl_t = 4


   integer(8), parameter :: BCS_DIRICHLET = 0  ! (E-field: electric wall)
   integer(8), parameter :: BCS_NEUMANN = 1    ! (E-field: magnetic wall)
   integer(8), parameter :: BCS_PERIODIC = 2

   integer(8), parameter :: FEM_FORMULATION_E = 1
   integer(8), parameter :: FEM_FORMULATION_H = 2


   integer(8), parameter :: NBERR_BAD_PERMUTATION     = -52
   integer(8), parameter :: NBERR_BAD_ADJACENCY       = -53
   integer(8), parameter :: NBERR_BAD_NODE_SEPARATION = -54




contains

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

        subroutine assert_or_die(pred, msg, ec)
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

            write(*,*) 'format', d_fmt
            write(buffer, d_fmt) val

            str = trim(buffer)
        end function int_2_str

        ! TODO: fix with just one call
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

            write(*,*) 'format', d_fmt
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



    end module numbatmod
