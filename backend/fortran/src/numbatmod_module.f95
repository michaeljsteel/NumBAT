

#include "numbat_decl.h"

module numbatmod

   implicit none

   integer, parameter :: STRINGLEN = D_STRINGLEN  ! defined in numbat_decl.h

   integer, parameter :: MAX_N_PTS = 250000
   integer, parameter :: MAX_N_ELTS = 120000

   integer, parameter :: MAX_LONG_ADJ = 2500000

   double precision, parameter :: D_PI = 3.141592653589793d0
   double precision, parameter :: D_ONE = 1.0d0
   double precision, parameter :: D_ZERO = 0.0d0
   complex*16, parameter :: C_IM_ONE = cmplx(0.0d0, 1.0d0, 8)

   integer*8, parameter :: nnodes_0 = 6
   integer*8, parameter :: nddl_t = 4


   integer*8, parameter :: BCS_DIRICHLET = 0  ! (E-field: electric wall)
   integer*8, parameter :: BCS_NEUMANN = 1    ! (E-field: magnetic wall)
   integer*8, parameter :: BCS_PERIODIC = 2

   integer*8, parameter :: FEM_FORMULATION_E = 1
   integer*8, parameter :: FEM_FORMULATION_H = 2




contains

   subroutine assert_no_larger_than(val, limit, location, msg, failco, errco, emsg)

      implicit none

      integer errco
      character emsg*(STRINGLEN)
      character location*(*), msg*(*)
      integer val, limit, failco

      if (val .ge. limit) then
         write(emsg,*) 'Failed limit check at ', location, '.  ', &
            'Expected ', msg, ',  but found values', val, limit
         errco = failco
      endif

      return
   end subroutine

end module numbatmod
