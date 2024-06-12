

#include "numbat_decl.h"

      module numbatmod

      implicit none

         integer, parameter :: STRINGLEN = D_STRINGLEN

         integer, parameter :: MAX_N_PTS = 250000
         integer, parameter :: MAX_N_ELTS = 120000

         integer, parameter :: MAX_LONG_ADJ = 2500000


      contains

         subroutine assert_no_larger_than(val, limit, location,
     *    msg, failco, errco, emsg)

            implicit none

            integer errco
            character emsg*(STRINGLEN)
            character location*(*), msg*(*)
            integer val, limit, failco

            if (val .ge. limit) then
               write(emsg,*) 'Failed limit check at ', location, '.  ',
     *         'Expected ', msg, ',  but found values', val, limit
               errco = failco
            endif

            return
         end subroutine

      end module numbatmod
